"""
AbstractQueue is an abstract type.
    There are three types:
    MyStack(): A Last In First Out Queue.
    FIFOQueue(): A First In First Out Queue.
    MyPriorityQueue(f,ext): Queue where items are sorted by f, (default <).
Each type supports the following methods and functions:
    append!(q,item)  -- add an item to the queue
    extend!(q,items) -- equivalent to: for item in items: append(q,item)
    pop!(q)          -- return the top item from the queue
    length(q)        -- number of items in q
"""
abstract type AbstractQueue{T} end

function extend!(Q::AbstractQueue, items)
    for item in items
        append!(Q, item)
    end
end

"""
Return an empty list, suitable as a Last-In-First-Out Queue.
"""
struct MyStack{T}
    e::DataStructures.Stack{T}
end
MyStack{T}() where {T} = MyStack{T}(DataStructures.Stack{T}())
append!(s::MyStack{T}, item::T) where {T} = push!(s.e, item)
Base.pop!(s::MyStack) = pop!(s.e)
Base.isempty(s::MyStack) = isempty(s.e)
Base.length(s::MyStack) = length(s.e)
Base.first(s::MyStack) = first(s.e)
Base.last(s::MyStack) = last(s.e)
Base.empty!(s::MyStack) = (empty!(s.e); s)
Base.:(==)(x::MyStack, y::MyStack) = x.e == y.e
function extend!(Q::MyStack, items)
    for item in items
        append!(Q, item)
    end
end

"""
A First-In-First-Out Queue.
"""
struct FIFOQueue{T}
    e::DataStructures.Queue{T}
end
FIFOQueue{T}() where {T} = FIFOQueue{T}(DataStructures.Queue{T}())
append!(q::FIFOQueue{T}, item::T) where {T} = DataStructures.enqueue!(q.e, item)
Base.pop!(q::FIFOQueue) = DataStructures.dequeue!(q.e)
Base.isempty(q::FIFOQueue) = isempty(q.e)
Base.length(q::FIFOQueue) = length(q.e)
Base.first(q::FIFOQueue) = first(q.e)
Base.last(q::FIFOQueue) = last(q.e)
Base.empty!(q::FIFOQueue) = (empty!(q.e); q)
Base.:(==)(x::FIFOQueue, y::FIFOQueue) = x.e == y.e
function extend!(Q::FIFOQueue, items)
    for item in items
        append!(Q, item)
    end
end

"""
A queue in which the minimum (or maximum) element (as determined by f)
is returned first. Keys of type T and priorities of type V.
"""
struct MyPriorityQueue{T, V}
    e::DataStructures.PriorityQueue{T, V}
    f::Any
    ext::Any
end

#check if it faster witn BinaryMinHeap or with a self-made implementation
MyPriorityQueue{T, V}(f, ext) where {T, V} =
    MyPriorityQueue{T, V}(DataStructures.PriorityQueue{T, V}(), f, ext)
append!(pq::MyPriorityQueue{T, V}, item::T) where {T, V} =
    DataStructures.enqueue!(pq.e, item, pq.f(item, pq.ext))
dequeue_pair!(pq::MyPriorityQueue) = DataStructures.dequeue_pair!(pq.e)
Base.pop!(pq::MyPriorityQueue) = DataStructures.dequeue!(pq.e)
Base.isempty(pq::MyPriorityQueue) = isempty(pq.e)
Base.length(pq::MyPriorityQueue) = length(pq.e)
Base.first(pq::MyPriorityQueue) = first(pq.e)
Base.last(pq::MyPriorityQueue) = last(pq.e)
Base.empty!(pq::MyPriorityQueue) = (empty!(pq.e); pq)
Base.:(==)(x::MyPriorityQueue, y::MyPriorityQueue) = x.e == y.e
function extend!(pq::MyPriorityQueue, items)
    for item in items
        append!(pq, item)
    end
end
