module DummyModule

export MyType, f

type MyType
    a::Int
end

f(x) = x^2+1

data = [ f(i) for i in 1:10 ]

println("DummyModule loaded on processor ",myid())

end

