// Robot/automaton product -> Spot parity game -> strategy and finite trace.
// Spot is an external GPLv3+ dependency, dynamically linked by this bridge.
#include <spot/parseaut/public.hh>
#include <spot/twaalgos/game.hh>
#include <spot/twaalgos/parity.hh>
#include <spot/twaalgos/isdet.hh>
#include <spot/twaalgos/hoa.hh>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <cstdint>
#include <limits>
#include <algorithm>

using U=uint32_t;
using V=uint64_t;
static U read32(std::istream& in) {
  unsigned char b[4]; if(!in.read(reinterpret_cast<char*>(b),4)) throw std::runtime_error("Truncated plant binary");
  return U(b[0])|(U(b[1])<<8)|(U(b[2])<<16)|(U(b[3])<<24);
}
static V read64(std::istream& in) { V a=read32(in),b=read32(in); return a|(b<<32); }
static void write32(std::ostream& out,U a) { for(int i=0;i<4;++i) out.put(char((a>>(8*i))&255)); }
static void write64(std::ostream& out,V a) { write32(out,U(a));write32(out,U(a>>32)); }
struct Plant {
  U nx,nu;
  std::vector<U> starts;
  std::vector<std::string> aps;
  std::vector<V> labels;
  std::vector<std::vector<U>> post;
  explicit Plant(const std::string& path) {
    std::ifstream in(path,std::ios::binary); char magic[8];in.read(magic,8);
    if(std::string(magic,8)!="ABCDTS01") throw std::runtime_error("Invalid plant binary magic");
    nx=read32(in);nu=read32(in);U ni=read32(in),na=read32(in);
    if(!nx || !nu || !ni || na>63) throw std::runtime_error("Invalid plant dimensions");
    for(U i=0;i<na;++i) { U len=read32(in);std::string name(len,' ');in.read(name.data(),len);aps.push_back(name); }
    for(U i=0;i<ni;++i) {U x=read32(in);if(x>=nx) throw std::runtime_error("Invalid initial state");starts.push_back(x);}
    for(U x=0;x<nx;++x) labels.push_back(read64(in));
    post.resize(V(nx)*nu);
    for(auto& successors:post) {U n=read32(in);successors.reserve(n);for(U i=0;i<n;++i) {U x=read32(in);if(x>=nx) throw std::runtime_error("Invalid successor");successors.push_back(x);}}
    if(!in || in.peek()!=EOF) throw std::runtime_error("Malformed plant binary");
  }
};
struct Transition {U q,priority;};
struct Node {U x,q;int action;};
struct Arc {U to,priority;};

// Independent universal verifier: an adversarial rejecting run exists iff
// some reachable cycle has even maximum priority. For every even c, search
// SCCs of the <=c subgraph for a c-edge contained in that SCC. Merely checking
// the maximum priority of whole SCCs is unsound (they can hide rejecting cycles).
static bool verify(const std::vector<std::vector<Arc>>& graph,U start,U colors) {
  U n=graph.size();std::vector<char> reachable(n);std::vector<U> todo{start};reachable[start]=1;
  for(size_t i=0;i<todo.size();++i) {
    U v=todo[i];if(graph[v].empty()) return false;
    for(auto e:graph[v]) if(!reachable[e.to]) {reachable[e.to]=1;todo.push_back(e.to);}
  }
  for(U c=0;c<colors;c+=2) {
    std::vector<char> seen(n);std::vector<U> order;
    for(U seed:todo) if(!seen[seed]) {
      std::vector<std::pair<U,size_t>> stack{{seed,0}};seen[seed]=1;
      while(!stack.empty()) {
        auto& top=stack.back();U v=top.first;
        while(top.second<graph[v].size() && (graph[v][top.second].priority>c || !reachable[graph[v][top.second].to])) ++top.second;
        if(top.second==graph[v].size()) {order.push_back(v);stack.pop_back();continue;}
        U w=graph[v][top.second++].to;
        if(!seen[w]) {seen[w]=1;stack.push_back({w,0});}
      }
    }
    std::vector<std::vector<U>> reverse(n);
    for(U v:todo) for(auto e:graph[v]) if(e.priority<=c && reachable[e.to]) reverse[e.to].push_back(v);
    std::vector<int> component(n,-1);int next=0;
    for(auto it=order.rbegin();it!=order.rend();++it) if(component[*it]<0) {
      std::vector<U> stack{*it};component[*it]=next;
      while(!stack.empty()) {U v=stack.back();stack.pop_back();for(U w:reverse[v]) if(component[w]<0) {component[w]=next;stack.push_back(w);}}
      ++next;
    }
    for(U v:todo) for(auto e:graph[v]) if(e.priority==c && component[v]==component[e.to]) return false;
  }
  return true;
}

int main(int argc,char** argv) {
  try {
    if(argc==2 && std::string(argv[1])=="--self-test") {
      if(verify({{{0,0},{1,1}},{{0,1}}},0,2)) throw std::runtime_error("Verifier missed a rejecting subcycle");
      if(!verify({{{1,0}},{{0,1}}},0,2)) throw std::runtime_error("Verifier rejected an accepting cycle");
      if(verify({{{1,3}},{{1,2}}},0,4)) throw std::runtime_error("Verifier counted a transient odd priority");
      if(!verify({{{0,1}},{{1,0}}},0,2)) throw std::runtime_error("Verifier counted unreachable rejection");
      if(verify({{}},0,2)) throw std::runtime_error("Verifier accepted a deadlock");
      std::cout<<"Five universal parity verifier cases passed\n";return 0;
    }
    if(argc!=8) throw std::runtime_error("Usage: bridge PLANT DPA_HOA RESULT_JSON POLICY_BIN HORIZON MAX_NODES MAX_EDGES");
    Plant plant(argv[1]);U horizon=std::stoul(argv[5]),max_nodes=std::stoul(argv[6]);V max_edges=std::stoull(argv[7]);
    auto dict=spot::make_bdd_dict();auto parsed=spot::parse_aut(argv[2],dict);
    if(parsed->format_errors(std::cerr) || !parsed->aut) throw std::runtime_error("Cannot parse DPA HOA");
    auto aut=parsed->aut;
    if(!spot::is_deterministic(aut)) throw std::runtime_error("Parity specification must be deterministic");
    // Normalize all parity conventions and implicit/uncolored priorities.
    aut=spot::colorize_parity_here(aut,true);
    aut=spot::change_parity_here(aut,spot::parity_kind_max,spot::parity_style_odd);
    aut=spot::colorize_parity_here(aut,true);
    bool maximum,odd;
    if(!aut->acc().is_parity(maximum,odd,true) || (aut->num_sets()>1 && !maximum) || !odd) {
      spot::print_hoa(std::cerr,aut);throw std::runtime_error("DPA parity normalization failed");
    }
    std::vector<std::pair<int,U>> variables;
    for(auto ap:aut->ap()) {
      auto it=std::find(plant.aps.begin(),plant.aps.end(),ap.ap_name());
      if(it==plant.aps.end()) throw std::runtime_error("Specification proposition has no plant label");
      variables.push_back({aut->register_ap(ap),U(it-plant.aps.begin())});
    }
    std::vector<V> masks=plant.labels;std::sort(masks.begin(),masks.end());masks.erase(std::unique(masks.begin(),masks.end()),masks.end());
    std::unordered_map<V,U> mask_index;for(U i=0;i<masks.size();++i) mask_index[masks[i]]=i;
    std::vector<Transition> delta(V(aut->num_states())*masks.size());
    for(U q=0;q<aut->num_states();++q) for(U i=0;i<masks.size();++i) {
      bdd cube=bddtrue;
      for(auto [var,bit]:variables) cube &= ((masks[i]>>bit)&1)?bdd_ithvar(var):bdd_nithvar(var);
      U matches=0;
      for(auto& e:aut->out(q)) if((e.cond & cube)!=bddfalse) {
        U count=0,color=0;for(U c:e.acc.sets()) {color=c;++count;}
        if(count!=1) throw std::runtime_error("DPA is not colored");
        delta[V(q)*masks.size()+i]={e.dst,color};++matches;
      }
      if(matches!=1) throw std::runtime_error("DPA incomplete or nondeterministic on an observed label");
    }
    auto step=[&](U q,U x){return delta[V(q)*masks.size()+mask_index.at(plant.labels[x])];};
    auto game=spot::make_twa_graph(dict);
    game->set_acceptance(aut->num_sets(),spot::acc_cond::acc_code::parity_max_odd(aut->num_sets()));
    spot::region_t owners;
    std::vector<Node> nodes;
    std::vector<U> controls;
    std::unordered_map<V,U> lookup;
    V edges=0;
    auto node=[&](Node n,bool controller) {
      if(nodes.size()>=max_nodes) throw std::runtime_error("Parity product resource limit: maximum nodes exceeded");
      U id=game->new_state();nodes.push_back(n);owners.push_back(controller);return id;
    };
    auto edge=[&](U from,U to,U color) {
      if(++edges>max_edges) throw std::runtime_error("Parity product resource limit: maximum edges exceeded");
      game->new_edge(from,to,bddtrue,{color});
    };
    auto control=[&](U x,U q) {
      V key=V(q)*plant.nx+x;auto it=lookup.find(key);if(it!=lookup.end()) return it->second;
      U id=node({x,q,-1},true);lookup[key]=id;controls.push_back(id);return id;
    };
    // Environment chooses the initial state: all configured starts must win.
    U initial=node({0,0,-2},false);game->set_init_state(initial);
    std::vector<U> starts;
    for(U x:plant.starts) {U q=step(aut->get_init_state_number(),x).q;U id=control(x,q);starts.push_back(id);edge(initial,id,0);}
    U rejecting=node({0,0,-3},false);edge(rejecting,rejecting,0);
    for(size_t pos=0;pos<controls.size();++pos) {
      U id=controls[pos],x=nodes[id].x,q=nodes[id].q;bool enabled=false;
      for(U u=0;u<plant.nu;++u) {
        auto& successors=plant.post[V(x)*plant.nu+u];if(successors.empty()) continue;
        enabled=true;U environment=node({x,q,int(u)},false);edge(id,environment,0);
        for(U xp:successors) {auto tr=step(q,xp);edge(environment,control(xp,tr.q),tr.priority);}
      }
      if(!enabled) edge(id,rejecting,0);
    }
    spot::set_state_players(game,std::move(owners));
    // Solve globally, then project only label-consistent automaton INITIAL
    // states. A nonempty raw product region is not sufficient for issuance.
    spot::solve_parity_game(game,true);
    auto& winners=spot::get_state_winners(game);auto& strategy=spot::get_strategy(game);
    std::vector<U> winning_starts;
    for(U id:starts) if(winners[id]) winning_starts.push_back(id);
    bool winning=!winning_starts.empty();
    U winning_products=0;for(U id:controls) if(winners[id]) ++winning_products;
    std::vector<int> action(nodes.size(),-1);
    std::vector<std::vector<Arc>> closed;
    bool verified=false;
    if(winning) {
      closed.resize(nodes.size());
      for(U v=0;v<nodes.size();++v) {
        bool controller=spot::get_state_player(game,v);
        if(controller && winners[v]) {
          U e=strategy[v];if(!e) throw std::runtime_error("Missing controller strategy on winning state");
          auto chosen=game->edge_storage(e);
          if(chosen.src!=v || nodes[chosen.dst].action<0) throw std::runtime_error("Invalid strategy/action mapping");
          action[v]=nodes[chosen.dst].action;
          U c=*chosen.acc.sets().begin();closed[v].push_back({chosen.dst,c});
        } else if(!controller) {
          for(auto& e:game->out(v)) closed[v].push_back({e.dst,*e.acc.sets().begin()});
        }
      }
      // Universally verify the strategy from EVERY projected winning start,
      // excluding losing candidate starts rather than requiring them to win.
      closed[initial].clear();
      for(U id:winning_starts) closed[initial].push_back({id,0});
      verified=verify(closed,initial,aut->num_sets());
      if(!verified) throw std::runtime_error("Independent universal parity verification failed");
      std::ofstream policy(argv[4],std::ios::binary);policy.write("ABCDPG01",8);
      write32(policy,aut->num_states());write32(policy,aut->num_sets());write32(policy,masks.size());write32(policy,controls.size());
      for(U id:controls) {write32(policy,nodes[id].x);write32(policy,nodes[id].q);write32(policy,U(action[id]+1));}
      for(V mask:masks) write64(policy,mask);
      for(auto tr:delta) {write32(policy,tr.q);write32(policy,tr.priority);}
      write32(policy,winning_starts.size());for(U id:winning_starts) {write32(policy,nodes[id].x);write32(policy,nodes[id].q);}
      if(!policy) throw std::runtime_error("Cannot export parity policy");
    }
    std::ofstream out(argv[3]);
    out<<"{\"winning\":"<<(winning?"true":"false")<<",\"verified\":"<<(verified?"true":"false")
       <<",\"game_nodes\":"<<nodes.size()<<",\"game_edges\":"<<edges<<",\"automaton_states\":"<<aut->num_states()
       <<",\"priorities\":"<<aut->num_sets()<<",\"reachable_product_states\":"<<controls.size()
       <<",\"winning_product_states\":"<<winning_products
       <<",\"candidate_initial_state_count\":"<<starts.size()
       <<",\"winning_initial_state_count\":"<<winning_starts.size()
       <<",\"winning_initial_states\":[";
    for(size_t i=0;i<winning_starts.size();++i) {if(i) out<<",";out<<nodes[winning_starts[i]].x+1;}
    out<<"]";
    if(winning) {
      U id=winning_starts.front();std::vector<U> xs{nodes[id].x+1},qs{nodes[id].q},us;
      for(U i=0;i<horizon;++i) {
        int u=action[id];if(u<0) throw std::runtime_error("Simulation reached undefined parity policy");
        U x=nodes[id].x,q=nodes[id].q;auto& successors=plant.post[V(x)*plant.nu+U(u)];
        if(successors.empty()) throw std::runtime_error("Simulation reached a deadlock");
        U xp=successors.front();for(U s:successors) if(s!=x) {xp=s;break;}
        U qp=step(q,xp).q;id=lookup.at(V(qp)*plant.nx+xp);
        xs.push_back(xp+1);qs.push_back(qp);us.push_back(U(u)+1);
      }
      auto array=[&](const char* name,const std::vector<U>& a) {out<<",\""<<name<<"\":[";for(size_t i=0;i<a.size();++i) {if(i) out<<",";out<<a[i];}out<<"]";};
      array("simulation_states",xs);array("simulation_memory",qs);array("simulation_inputs",us);
      out<<",\"simulation_steps\":"<<horizon<<",\"simulation_status\":\"horizon_reached\"";
    }
    out<<"}\n";if(!out) throw std::runtime_error("Cannot write parity result");
    return 0;
  } catch(const std::exception& e) {std::cerr<<e.what()<<"\n";return 2;}
}
