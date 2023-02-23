using Plots
using Primes
using LaTeXStrings
#import Pkg; Pkg.add("Plots"); Pkg.add("Primes")



# <Summary> 
# This Code generates to plot SL(2,Z) in ram.
# PlotGpElements
# PlotWeylsCriterium
#</Summary>

#This to stop making svg files.
default(fmt=:png)
###############################################################################
# ADD THIS BLOCK TO REMOVE SVG FROM LIST OF "DISPLAYABLE_MIMES":
pos = findfirst((x)->(x=="image/svg+xml"), VSCodeServer.DISPLAYABLE_MIMES)
if !isnothing(pos)
    popat!(VSCodeServer.DISPLAYABLE_MIMES, pos)
    println("Popped!")
end
###############################################################################

# <Summary>Permutions return a list of the isometric Permutations of a Copirme </Summary>
struct SLPoint
    a::Int
    b::Int
    c::Int
    d::Int
end

# <Summary>Permutions return a list of the isometric Permutations of a Copirme </Summary>
struct CoPrime
    m::Int
    n::Int
end

# <Summary>Give the Determinant of the SL Point </Summary>
# <param name="point"> SLPoint</param>
# <return>int</return>
function Determinant(point::SLPoint)
    return point.a*point.d - point.b*point.c 
 end

 # <Summary>Calculate the Fabenius Norm </Summary>
# <param name="point"> SLPoint </param>
# <return>int</return>
function Norm(point::SLPoint)
    return point.a^2 + point.b^2 + point.c^2 + point.d^2
end


# <Summary>Permutions return a list of the isometric Permutations of a Copirme </Summary>
# <param name="CoP">Coprime</param>
# <return>A list for Copirmes</return>
function Permutions(CoP::CoPrime)
    #Since (0,1) and (1,1) have some symetric we use unique we might just dont use it on those...    
    return unique([CoPrime(CoP.m,CoP.n),
    CoPrime(-CoP.m,CoP.n),
    CoPrime(CoP.m,-CoP.n),
    CoPrime(-CoP.m,-CoP.n),
    CoPrime(CoP.n,CoP.m),
    CoPrime(-CoP.n,CoP.m),
    CoPrime(CoP.n,-CoP.m),
    CoPrime(-CoP.n,-CoP.m)])
end

# <Summary>Recursive function that, append Postive Coprime Pairs to a vector</Summary>
# <param name="Co">Coprime</param>
# <param name="CoPrimeList">Vector of CoPrime</param>
# <param name="MAX">Int</param>
# <return></return>
function GeneratingPostiveCoprimePairs(Co::CoPrime,CoPrimeList::Vector{CoPrime},MAX)
    if  2*Co.m-Co.n <= sqrt(MAX - Co.m^2)
        push!(CoPrimeList,CoPrime(2*Co.m-Co.n,Co.m))
        GeneratingPostiveCoprimePairs(CoPrime(2*Co.m-Co.n,Co.m),CoPrimeList,MAX)
    end
    if 2*Co.m+Co.n  <= sqrt(MAX - Co.m^2)
        push!(CoPrimeList,CoPrime(2*Co.m+Co.n,Co.m))
        GeneratingPostiveCoprimePairs(CoPrime(2*Co.m+Co.n,Co.m),CoPrimeList,MAX)
    end
    if Co.m+2*Co.n  <= sqrt(MAX - Co.n^2)
        push!(CoPrimeList,CoPrime(Co.m+2*Co.n,Co.n))
        GeneratingPostiveCoprimePairs(CoPrime(Co.m+2*Co.n,Co.n),CoPrimeList,MAX)
    end 
end

# <Summary>
# Recursive function that, Generatin Postive Coprime Pairs (n,m)
# such that n^2+m^2 \leq MAX
# </Summary>
# <param name="MAX">Int</param>
# <param name="All"> if (1,1) and (1,0) is include</param>
# <return>Vector of Postive Coprimes pairs</return>
function GeneratingPostiveCoprimePairs(MAX, All = true)
    if All
        append!(CoPrimeList,[CoPrime(1,1)])
        append!(CoPrimeList,[CoPrime(1,0)])
    end
    
    CoPrimeList = Vector{CoPrime}()

    push!(CoPrimeList,CoPrime(2,1))
    GeneratingPostiveCoprimePairs(CoPrime(2,1),CoPrimeList,MAX)

    push!(CoPrimeList,CoPrime(3,1))
    GeneratingPostiveCoprimePairs(CoPrime(3,1),CoPrimeList,MAX)

    return CoPrimeList
end

#Get a range of Solutions for a Quadratic Equation
function QuadraticSolutions(a,b,c,d,M,m = 0)
    A = convert(Int128, a)
    B = convert(Int128, b)
    C = convert(Int128, c)
    D = convert(Int128, d)

    Mdet = (-A^4 - 2*A^2*B^2 + A^2*C^2 + 2*A*B*C*D - B^4 + B^2*D^2 + M*A^2 + M*B^2)
    if Mdet <= 0
        return 0:0
    end
    Mres1 = (A*C + B*D + sqrt(Mdet))/(A^2 + B^2)
    Mres2 = (A*C + B*D - sqrt(Mdet))/(A^2 + B^2)
    
    if 0 < m < M 
        mdet = (-A^4 - 2*A^2*B^2 + A^2*C^2 + 2*A*B*C*D - B^4 + B^2*D^2 + m*A^2 + m*B^2)
        if mdet <= 0
            return 0:0
        end
        mres1 = (A*C + B*D + sqrt(mdet))/(A^2 + B^2)
        mres2 = (A*C + B*D - sqrt(mdet))/(A^2 + B^2)
        return [floor(Int,min(Mres1, Mres2)): floor(Int,min(mres1, mres2));
                 ceil(Int,max(mres1, mres2)):  ceil(Int,max(Mres1, Mres2))]
    end

    return floor(Int,min(Mres1, Mres2)) : ceil(Int,max(Mres1, Mres2))
end

function BezoutGpElemtents(CoP::CoPrime, MAX, MIN=0,f::Function =_ -> true)
    return BezoutGpElemtents(CoP.m,CoP.n,MAX,MIN,f)
end


function BezoutGpElemtents(a,b,MAX,MIN=0,f::Function =_ -> true)
    SLPoints =  Vector{SLPoint}()
    gcd, d0, c0 = gcdx(a, -b)
            if gcd == 1
                c = c0
                d = d0
                for n in QuadraticSolutions(a,b,c0,d0,MAX, MIN)
                    d = d0 - n * b
                    c = c0 - n * a
                    NewSLPoint = SLPoint(a,b,c,d)
                    if Determinant(NewSLPoint) == 1 && MIN <=Norm(NewSLPoint) < MAX && f(NewSLPoint)
                        push!(SLPoints,NewSLPoint)
                    end
                end
            end
    return SLPoints
end

function BezoutGpElemtents(MAX, f::Function = _ -> true ,
                                g::Function = point::SLPoint -> point,
                                h::Function = CoP::CoPrime -> Permutions(CoP) )
    SLPoints =  Vector{SLPoint}()
    coPrimes = GeneratingPostiveCoprimePairs(MAX,false)

    for CoP in unique([h(CoPrime(1,1));h(CoPrime(0,1))])
        append!(SLPoints,g.(BezoutGpElemtents(CoP.m,CoP.n,MAX,0,f)))
    end

    for CoP in coPrimes
        for nCoP in h(CoP)
            append!(SLPoints,g.(BezoutGpElemtents(nCoP.m,nCoP.n,MAX,0,f)))
        end
    end
    return SLPoints
end

function WeylsCriterium(v,m=1)
    function NumberExp(base,input,m=1)
        if m == 0 
            return base + 1
        end
        return base+exp(input*1im*m)
    end
    
    vector = Vector{Number}()
    if m == 0 
        zero::Number = 0.0 
        vector = push!(vector,zero)
    end

    append!(vector,v)
    res = collect(Iterators.accumulate((y,x) -> NumberExp(y,x,m),vector))
    return res
end


function Poincare(z::Complex)
    return (z- 1im )/(z + 1im)
end

function Mobius(point::SLPoint)
    return (point.a*1im + point.b)/(point.c*1im + point.d)
end

function Inverse(point::SLPoint)
    return SLPoint(point.d,-point.b,-point.c,point.a)
end


function ToReal2(z::Complex)
    return (real(z),imag(z))
end

function ToRealHyb2(z::Complex)
    return (real(z)*(1-abs2(z))^(-1),imag(z)*(1-abs2(z))^(-1))
end



function ArgeumtOfComplex(point::SLPoint)
    z = Poincare(Mobius(point))
    iz = Poincare(Mobius(Inverse(point)))
    return (angle(z), angle(iz))
end



function PlotGpElements(MAX,header="",f::Function = _ -> true ,
                            g::Function = point::SLPoint -> point,
                            h::Function = CoP::CoPrime -> Permutions(CoP))
    GpElemts = BezoutGpElemtents(MAX,f,g,h)
    y1 = Mobius.(GpElemts)
    y2 = Poincare.(y1)

    y4 = (Mobius.(Inverse.(GpElemts)))
    y5 = Poincare.(y4)

    y7 = ArgeumtOfComplex.(GpElemts)

    println("Plotter")

    p1 = scatter(ToReal2.(y1), markersize=3,ma=0.3,label="", ylab="       i"
    , yguidefontrotation=-90,title="γ(i)")
    p2 = scatter(ToReal2.(y2), ylab="         i"
    , yguidefontrotation=-90,aspect_ratio = 1,markersize=3,ma=0.3,label="",title="Cayley Transformation")
    p3 = scatter(ToRealHyb2.(y2), aspect_ratio = 1,markersize=3,ma=0.3,label="",title="1/(1-r²)", yaxis=false,xaxis=false)
    
    p4 = scatter(ToReal2.(y4), ylab="       i"
    , yguidefontrotation=-90, markersize=3,ma=0.3,label="",title="γ⁻¹(i)")
    p5 = scatter(ToReal2.(y5), ylab="         i"
    , yguidefontrotation=-90,aspect_ratio = 1,markersize=3,ma=0.3,label="",title="Cayley Transformation")
    p6 = scatter(ToRealHyb2.(y5), aspect_ratio = 1,markersize=3,ma=0.3,label="",title="1/(1-r²)", yaxis=false,xaxis=false)
    
    p7 = scatter(y7, aspect_ratio = 1,markersize=3,ma=0.3,label="",xlabel="Arg(γ(i))",ylabel="\n Arg(γ⁻¹(i))",title="Cartan’s decomposition", xticks=([-pi,-pi/2 ,0, pi/2, pi], ["-π","-π/2","0","π/2", "π"]),yticks=([-pi,-pi/2 ,0, pi/2, pi], ["-π","-π/2","0","π/2", "π"]))
   
    lifted = replace("$(log10(MAX))","0"=>"⁰ ","1"=>"¹","2"=>"²","3"=>"³","4"=>"⁴","5"=>"⁵","6"=>"⁶","7"=>"⁷","8"=>"⁸","9"=>"⁹","."=>"ᐧ")
    Header = "γ ∈ "* header * ", ‖γ‖ ≤ 10"*lifted;

    l = @layout [ grid(2,3)
             b{0.4h} ]  


return  plot(p1, p2, p3, p4, p5, p6, p7, layout=l, legend=false,plot_title=Header
,sizes = (600,800))

end

# <Summary>SortSLPoints </Summary>
# <param name="GpElemtents">List of SL Points</param>
# <param name="MAX"> Vector of CoPrime</param>
# <param name="f">Sort Criterium</param>
# <return></return>
function SortSLPoints(GpElemtents, MAX, f::Function = Norm)
    #[[SLPoint(1,1,1,1),SLPoint(1,1,1,1)]]
    #Elements = sort(GpElemtents, by=f)
    
    List = Vector{Vector{SLPoint}}()
    for _ in range(1,MAX)
        push!(List,[])
    end
    for point in GpElemtents
        push!(List[Norm(point)], point)
    end

    return List
end

# <Summary>PlotWeylsCriterium </Summary>
# <param name="header">Coprime</param>
# <param name="plotHead">Vector of CoPrime</param>
# <param name="f">Filter - e.g. point -> isprime(Norm(point))</param>
# <param name="g">Int</param>
# <param name="h">Int</param>
# <param name="k">(x,X,N)  </param>
# <param name="l">(x,X,N)  </param>
# <return></return>
function PlotWeylsCriterium(MAX,M=0, 
    header="",
    plotHead="",
    xAxisLabel = "",
    yAxisLabel = "",
    f::Function = _ -> true ,
    g::Function = point::SLPoint -> point,
    h::Function = CoP::CoPrime -> Permutions(CoP),
    k::Function = X -> 6X,
    l::Function = (x , X, N) -> x
    )
    GpElemtents = BezoutGpElemtents(MAX,f,g,h)
    
    yLabel = yAxisLabel
    xLabel = xAxisLabel
    
    x1= range(0, MAX, length=100)
    y1 = @. k(x1)
    show = plot(x1,y1, 
        label=plotHead, 
        title=header,
        yaxis=yLabel, 
        xaxis=xLabel, 
        sizes = (1200,400)
        )
    sorted = SortSLPoints(GpElemtents,MAX, Norm)
    res::Number = 0.0;
    PlotsPoint = Vector{Float64}()
    NumberOfElements = 0;
    for SLPs in sorted
        if length(SLPs) != 0
            angs = angle.(Poincare.(Mobius.(SLPs)))
            res += last(WeylsCriterium(angs, M))
            NumberOfElements += length(SLPs)
        end
        push!(PlotsPoint, l(abs.(res),length(PlotsPoint)+1,NumberOfElements))    
    end
    plot!(PlotsPoint,label="")
    p1 = plot(ticks = false,axis=([], false))
    l = @layout [a{0.01w} grid(1,1)]  

return  plot(p1, show, layout=l,sizes = (1000,400))
end

#=
#Exampel
PlotWeylsCriterium(10^3,
    1, 
    "Weyl's Criterium for γ ∈ SL(2,ℤ) and ‖γ‖ is prime ",
    "",
    "X \n",
    L"|\dfrac{1}{P(X)} \sum_{‖γ‖_{F^2} ≤ X} e_1(arg(γ(i)))|",
    point -> isprime(Norm(point)) ,
    point::SLPoint -> point,
    CoP::CoPrime -> Permutions(CoP),
    X -> 0,
    (x , X, N) -> x/N
    )
                                                         =#
