import itertools
from pathlib import Path
import random
import sys
import unittest

from build_bridge import build
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"ParitySynthesis"))
from test_native import NativeTests as OriginalHarness


def fixed_policy_wins(labels,adjacency,start,actions,formula):
    reachable={start};todo=[start]
    for x in todo:
        successors=adjacency[x][actions[x]]
        if not successors:return False
        for xp in successors:
            if xp not in reachable:reachable.add(xp);todo.append(xp)
    if formula=="G blue":return all(labels[x] for x in reachable)
    def cycles_good(path):
        for xp in adjacency[path[-1]][actions[path[-1]]]:
            if xp in path:
                cycle=path[path.index(xp):]
                if formula=="F G blue" and not all(labels[x] for x in cycle):return False
                if formula=="G F blue" and not any(labels[x] for x in cycle):return False
            elif not cycles_good(path+[xp]):return False
        return True
    return all(cycles_good([x]) for x in reachable)


class RegionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.harness=OriginalHarness()
        cls.harness.bridge=build()

    def test_every_projected_state_against_exhaustive_strategy_oracle(self):
        rng=random.Random(212026)
        comparisons=0
        for _ in range(120):
            n=rng.randint(1,4)
            labels=[rng.randrange(2) for _ in range(n)]
            adjacency=[[[xp for xp in range(n) if rng.random()<.48] for _ in range(2)] for _ in range(n)]
            for formula in ("G blue","G F blue","F G blue"):
                expected={x+1 for x in range(n) if any(fixed_policy_wins(labels,adjacency,x,a,formula)
                    for a in itertools.product(range(2),repeat=n))}
                result=self.harness.solve(labels,adjacency,list(range(n)),formula)
                self.assertEqual(set(result["winning_initial_states"]),expected,(labels,adjacency,formula))
                self.assertEqual(result["winning"],bool(expected))
                self.assertEqual(result["winning_initial_state_count"],len(expected))
                comparisons+=1
        print("Exact projected-region comparisons:",comparisons)

    def test_some_initial_states_win_and_others_lose(self):
        result=self.harness.solve([1,0],[[[0]],[[1]]],[0,1],"F G blue")
        self.assertTrue(result["winning"] and result["verified"])
        self.assertEqual(result["winning_initial_states"],[1])
        self.assertEqual(result["candidate_initial_state_count"],2)


if __name__=="__main__":
    # Avoid unittest collecting the imported legacy harness, whose semantics
    # deliberately require all configured starts rather than any winning start.
    unittest.TextTestRunner().run(unittest.defaultTestLoader.loadTestsFromTestCase(RegionTests)).wasSuccessful() or sys.exit(1)
