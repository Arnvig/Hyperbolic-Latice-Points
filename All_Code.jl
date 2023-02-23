using Plots
using Primes
using LaTeXStrings

using BenchmarkTools


default(fmt=:png)
#import Pkg; Pkg.add("Plots"); Pkg.add("Primes")
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


# <Summary>Recursive function that, that  </Summary>
# <param name="Co">Coprime</param>
# <param name="CoPrimeList">Vector of CoPrime</param>
# <param name="Max">Int</param>
# <return></return>
function GeneratingPostiveCoprimePairs(Co::CoPrime,CoPrimeList::Vector{CoPrime},Max)
    if  2*Co.m-Co.n <= sqrt(Max - Co.m^2)
        push!(CoPrimeList,CoPrime(2*Co.m-Co.n,Co.m))
        GeneratingPostiveCoprimePairs(CoPrime(2*Co.m-Co.n,Co.m),CoPrimeList,Max)
    end
    if 2*Co.m+Co.n  <= sqrt(Max - Co.m^2)
        push!(CoPrimeList,CoPrime(2*Co.m+Co.n,Co.m))
        GeneratingPostiveCoprimePairs(CoPrime(2*Co.m+Co.n,Co.m),CoPrimeList,Max)
    end
    if Co.m+2*Co.n  <= sqrt(Max - Co.n^2)
        push!(CoPrimeList,CoPrime(Co.m+2*Co.n,Co.n))
        GeneratingPostiveCoprimePairs(CoPrime(Co.m+2*Co.n,Co.n),CoPrimeList,Max)
    end 
end

#Generates All CoPrime(m,n) with m>n 
function GeneratingPostiveCoprimePairs(Max, All = true)
    if All
        append!(CoPrimeList,[CoPrime(1,1)])
        append!(CoPrimeList,[CoPrime(1,0)])
    end
    
    CoPrimeList = Vector{CoPrime}()

    push!(CoPrimeList,CoPrime(2,1))
    GeneratingPostiveCoprimePairs(CoPrime(2,1),CoPrimeList,Max)

    push!(CoPrimeList,CoPrime(3,1))
    GeneratingPostiveCoprimePairs(CoPrime(3,1),CoPrimeList,Max)

    return CoPrimeList
end


function ToFold(CoP::CoPrime, MAX)
    function InternalFoldFunction(x)
        return BezoutGpElemtents(x, MAX)
    end
    return reduce(vcat,InternalFoldFunction.(Permutions(CoP)))
end

function FoldPostiveCoprimePairs(MAX, f :: Function = ToFold , All = true)
    SLList = Vector{SLPoint}()
    append!(SLList,f(CoPrime(2,1),MAX))
    FoldPostiveCoprimePairs(CoPrime(2,1),SLList,f,MAX)
    append!(SLList,f(CoPrime(3,1),MAX))
    FoldPostiveCoprimePairs(CoPrime(3,1),SLList,f,MAX)

    if All
        append!(SLList,f(CoPrime(1,1),MAX))
        append!(SLList,f(CoPrime(1,0),MAX))
    end
    
    return SLList
end

function FoldPostiveCoprimePairs(Co::CoPrime,SLList::Vector{SLPoint}, f::Function,MAX)
    if  2*Co.m-Co.n <= sqrt(MAX - Co.m^2)
        append!(SLList,f(CoPrime(2*Co.m-Co.n,Co.m),MAX))
        FoldPostiveCoprimePairs(
            CoPrime(2*Co.m-Co.n,Co.m),SLList,f,MAX)
    end
    if 2*Co.m+Co.n  <= sqrt(MAX - Co.m^2)
        append!(SLList,f(CoPrime(2*Co.m+Co.n,Co.m),MAX))
        FoldPostiveCoprimePairs(CoPrime(2*Co.m+Co.n,Co.m),SLList,f,MAX)
    end
    if Co.m+2*Co.n  <= sqrt(MAX - Co.n^2)
        append!(SLList,f(CoPrime(Co.m+2*Co.n,Co.n),MAX))
        FoldPostiveCoprimePairs(CoPrime(Co.m+2*Co.n,Co.n),SLList,f,MAX)
    end 
end


function BezoutGpElemtents(CoP::CoPrime, MAX, MIN=0,f::Function =_ -> true)
    return BezoutGpElemtents(CoP.m,CoP.n,MAX,MIN,f)
end

function Bezout(CoP::CoPrime)
    Bezout(CoP.n,CoP.m)
end

function Bezout(a,b)
    gcd, d0, c0 = gcdx(a, -b)
    SLPoints =  Vector{SLPoint}()
    gcd, d0, c0 = gcdx(a, -b)
        if gcd == 1
            k = (a*c0 + b*d0)/(a^2 + b^2)
            d1 = d0 - floor(k) * b
            c1 = c0 - floor(k) * a
            d2 = d0 - ceil(k) * b
            c2 = c0 - ceil(k) * a

            if c1^2 + d1^2 <= c2^2 + d1^2 
                push!(SLPoints,SLPoint(a,b,c1,d1))
            else 
                push!(SLPoints,SLPoint(a,b,c2,d2))
            end 
        end
    return SLPoints
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
    #println(res)
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

function Transpose(point::SLPoint)
    return SLPoint(point.d,point.b,point.c,point.a)
end

function DiagonalTranspose(point::SLPoint)
    return SLPoint(point.a,point.c,point.b,point.d)
end

function ToReal2(z::Complex)
    return (real(z),imag(z))
end

function ToRealHyb2(z::Complex)
    return (real(z)*(1-abs2(z))^(-1),imag(z)*(1-abs2(z))^(-1))
end

function ToTuple(Co::CoPrime)
    return (Co.m,Co.n)
end

function ArgeumtOfComplex(point::SLPoint)
    z = Poincare(Mobius(point))
    iz = Poincare(Mobius(Inverse(point)))
    return (angle(z), angle(iz))
end

function TransformedArgeumtOfComplex(point::SLPoint)
    z = Poincare(Mobius(point))
    iz = Poincare(Mobius(Inverse(point)))
    return (angle(iz), angle(z))
end

#Since we do not want to get error we do it a fun way.
function issquare(x)
    abs(x)==floor(sqrt(abs(x)))^2 
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


function CoPrimesEliminator(coPrime::CoPrime)
    if coPrime == CoPrime(1,1) || coPrime == CoPrime(0,1)
        return []
    end
    return [coPrime]
end



function TestPermutions(CoP::CoPrime)
    if CoP == CoPrime(1,1) || CoP == CoPrime(0,1)
        return []
    end
    return [CoPrime(CoP.m,CoP.n),CoPrime(CoP.n,CoP.m)]
end







function Plot_WeylsCriterium(max, M =10 ,f::Function = _::SLPoint -> true)
    x = angle.(Poincare.(Mobius.(BezoutGpElemtents(max,f))))
    y = abs.(WeylsCriterium(x))
    x1= range(0, length(x), length=100)
    y1 = @. (6/(pi^2))*x1
    #
    show = plot(x1,y1, label="x/Log(x)")
    println("starting ploting")
    for m in range(0,M) 
        y = abs.(WeylsCriterium(x,m))
        plot!(y,label="",
                color=cgrad(:thermal, rev = true)[floor(Int,255*m/M)+1])
        println(m)
    end
    return show
end

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
# <param name="k">Int</param>
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
    println("Last Point")
    print(PlotsPoint[lastindex(PlotsPoint)])
    p1 = plot(ticks = false,axis=([], false))
    l = @layout [a{0.01w} grid(1,1)]  


return  plot(p1, show, layout=l,sizes = (1000,400))

    #return show
end

# <Summary>Permutions return a list of the isometric Permutations of a Copirme </Summary>
# <param name="Max">int </param>
# <param name="M">int What power, default 1, 0 counts points  </param>
# <param name="f"> Filter for SLPoint </param>
# <param name="h"> Filter for SLPoint </param>
# <return>Plot</return>

function PlotWeylsCriterium(MAX, M =1 ,
    f::Function = _::SLPoint -> true, 
    h::Function = x -> 6*x, 
    g::Function = (x,y) -> x )
    #x = angle.(Poincare.(Mobius.(BezoutGpElemtents(MAX,f))))
    #y = abs.(WeylsCriterium(x))
    x1= range(0, MAX, length=100)
    y1 = @. h(x1)
    #x1,y1, label="x/Log(x)"
    show = plot(x1,y1, label="x/Log(x)")
    println("starting ploting")
    res = 0.0;
    PlotsPoint = Vector{Number}()

    GpElemtents = BezoutGpElemtents(MAX,f)
    sorted = SortSLPoints(GpElemtents,MAX, Norm)
    
    for SLPs in sorted
        if length(SLPs) != 0
            angs = angle.(Poincare.(Mobius.(SLPs)))
            res += last(WeylsCriterium(angs, M))
        end
        push!(PlotsPoint, g(abs.(res),length(PlotsPoint)+1))    
    end
    plot!(PlotsPoint,label="")
    print(PlotsPoint[lastindex(PlotsPoint)])
    return show
end



#PlotWeylsCriterium(10^4,0, SLPoint -> 0<SLPoint.b<SLPoint.a, x -> 2/3*x)
#PlotWeylsCriterium(10^4,0, SLPoint -> 0<SLPoint.a<SLPoint.b, x -> 2/3*x)
#PlotWeylsCriterium(10^4,0, SLPoint -> SLPoint.a<SLPoint.b<0, x -> 2/3*x)
#PlotWeylsCriterium(10^4,0, SLPoint -> SLPoint.b<SLPoint.a<0, x -> 2/3*x)

#PlotWeylsCriterium(10^2,0, SLPoint -> SLPoint.b==SLPoint.a, x -> 0)

#PlotGpElements(10^4,"GCD(N,M)=1 and M<N",_ -> true, point -> point, CoPrimesEliminator)

function NotCoPrimesEliminator(coPrime::CoPrime)
    if coPrime == CoPrime(1,1) 
        return [CoPrime(1,1)]
    elseif coPrime == CoPrime(0,1)
        return [CoPrime(0,1)]
    end
    return []
end

function TimeGlasPrimesEliminator(coPrime::CoPrime)
    if coPrime == CoPrime(1,1) 
        return [CoPrime(1,1)]
    elseif coPrime == CoPrime(0,1)
        return []
    end
    return []
end

function FireKuglerCoPrimesEliminator(coPrime::CoPrime)
    if coPrime == CoPrime(1,1) 
        return [CoPrime(1,1), CoPrime(1,-1)]
    elseif coPrime == CoPrime(0,1)
        return [CoPrime(0,1), CoPrime(1,0)]
    end
    return []
end

function CoPrimesPermutions(coPrime::CoPrime)
    if coPrime == CoPrime(1,1) 
        return Permutions(coPrime)
    elseif coPrime == CoPrime(0,1)
        return Permutions(coPrime)
    end
    return []
end


#PlotWeylsCriterium(10^3,0,"", _ -> true ,    point-> point,CoPrimesEliminator, x -> (2/3)*x)

function StoreO(x)
    if x < 2
        return 0
    elseif x < 5
        2*sqrt(x-2)
    end
    return 2*sqrt(x-2)+2*sqrt(x-5)

end

function TestElements(MAX)
    len = length(BezoutGpElemtents(MAX, _ -> true, point -> point, NotCoPrimesEliminator))
    isIt = StoreO(MAX)

    return isIt - len
end

#PlotWeylsCriterium(10^6,0,"(1,1) & (0,1)", "2*sqrt(x-2)+2*sqrt(x-5)", _ -> true , point-> point ,NotCoPrimesEliminator,StoreO)
#PlotGpElements(10^4,"GCD(N,M)=1 and M<N",_ -> true, point -> point, NotCoPrimesEliminator)

#TestElements(10^7)

#TestElements(10^5)

function CoPrimesEliminatorAndPermutions(coPrime::CoPrime)
    if coPrime == CoPrime(1,1) || coPrime == CoPrime(0,1)
        return []
    end
    return Permutions(coPrime)
end




#PlotWeylsCriterium(10^7,0,"SL / (1,1,k,k+1) and (0,1,-1,k) there's Permutations", "f(x)=6x", _ -> true , point-> point ,CoPrimesEliminatorAndPermutions, x -> 6*x)



#PlotWeylsCriterium(10^3,0,"Permutions((1,1)) and Permutions((0,1))", "8*sqrt(x)", _ -> true , point-> point ,CoPrimesPermutions,x -> 8*sqrt(x))

function PlotCoPrime(MAX)
    points = GeneratingPostiveCoprimePairs(MAX,false)
    append!(points, [CoPrime(0,1),CoPrime(1,1)])
    points = Permutions.(points)
    points= reduce(vcat,points)
    
     

    function ToXYZ(CoP::CoPrime)
        return (CoP.m,CoP.n,length(BezoutGpElemtents(CoP,MAX)))
    end
    #Permutions

    points = ToXYZ.(points)
    points = filter(point -> 0<point[3]<sqrt(sqrt(MAX)), points)

    y = getindex.(points,1)
    x = getindex.(points,2)
    z = getindex.(points,3)

    scatter(x,y, title="CoPrimes", 
        label="", 
        yaxis="n", 
        xaxis="m", 
        zlabel="k",
        zcolor=-z
        ,sizes = (1200,1000))
end

function Plot3DCoPrime(MAX)
    points = GeneratingPostiveCoprimePairs(MAX,false)
    append!(points, [CoPrime(0,1),CoPrime(1,1)])
    points = Permutions.(points)
    points= reduce(vcat,points)
    #GeneratingPostiveCoprimePairs(MAX,false)
    function ToX(CoP::CoPrime)
        return CoP.m
    end
    function ToY(CoP::CoPrime)
        return CoP.n
    end
    function ToZ(CoP::CoPrime)
        return length(BezoutGpElemtents(CoP,MAX))
    end
    function ToZerror(CoP::CoPrime)
        return (length(BezoutGpElemtents(CoP,MAX)),0)
    end
    #Permutions

    y = ToY.(points)
    x = ToX.(points)
    z = ToZ.(points)
    zerr = ToZerror.(points)

    scatter(x,y,z, title="CoPrimes", 
        label="", 
        yaxis="n", 
        xaxis="m", 
        zlabel="k",
        zerror=zerr,
        zcolor=-z,
        zlims=(1,maximum(z))
        ,sizes = (1200,1000))
end



function plotGPSolutions(MAX)

    points = GeneratingPostiveCoprimePairs(MAX,false)
    append!(points, [CoPrime(0,1),CoPrime(1,1)])

    #x1,y1, label="x/Log(x)"

    show = scatter()
    plot!( label="")
    println("starting ploting")
    PlotsPoint = Vector{Float64}()

    function Sorter(CoP::CoPrime)
        return CoP.m^2+CoP.n^2
    end
    sorted = sort(points, by=Sorter)
    
    for SLPs in sorted
        push!(PlotsPoint, length(BezoutGpElemtents(SLPs,MAX)))    
    end
    scatter!(PlotsPoint,label="")
    return plot(show)
end
#PlotCoPrime(10^3.5)

#Plot3DCoPrime(10^3)

#PlotCoPrime(10^5)

#plotGPSolutions(10)



#PlotWeylsCriterium(10^1,0,"SL", "6*x", _ -> true , point-> point ,Permutions,x -> 6*x)

#PlotWeylsCriterium(10^3,0,"SL", "6*x", _ -> true , point-> point ,Permutions,x -> 6*x )

function test10(x,y)
    return (x-6*(y))/((y)^(0.5))
end
#(x,y)-> (x-6*y)/(y^(0.5))

#PlotWeylsCriterium(10^3, 0, _ -> true, x->0, test10)

#PlotWeylsCriterium(10^6, 0, point -> isprime(Norm(point)), x->0, (x,y)-> (x - 4.3*y/log(y)))
#PlotWeylsCriterium(10^6, 0, point -> isprime(Norm(point)-2), x->0, (x,y)-> (x-0.831*y))


#PlotWeylsCriterium(10^6, 1, point -> isprime(Norm(point)-2), x->0, (x,y)-> (x/y))

#PlotWeylsCriterium(10^6, 1, point -> isprime(Norm(point)), x->0, (x,y)-> (x/y))

#=
PlotWeylsCriterium(10^3,
    1, 
    "",
    "",
    "",
    "",
    point -> isprime(Norm(point)-2) ,
    point::SLPoint -> point,
    CoP::CoPrime -> Permutions(CoP),
    X -> 1,
    (x , X, N) -> x/N
    )
=# 

#=PlotWeylsCriterium(10^7
,0,
"SL / (1,1,k,k+1) and (0,1,-1,k) there's Permutations", 
"f(x)=6x", _ -> true ,
 point-> point ,
 CoPrimesEliminatorAndPermutions, x -> 6*x)
 
=#
 
#PlotGpElements(10^5,"0<a<b ", point -> 0<point.a<point.b , point -> point)
#PlotGpElements(10^3.5,"0<b<a", point -> 0<point.b<point.a , point -> point)
#PlotGpElements(10^3.5,"0<b<-a", point -> 0<point.b<-point.a , point -> point)
#PlotGpElements(10^3.5,"0<a<-b ", point -> 0 <point.a < -point.b , point -> point)


#=
PlotGpElements(10^3.5,"0<a<-b ", point -> 0 <point.a < -point.b , point -> point)
PlotWeylsCriterium(10^5,
    1, 
    "",
    "",
    "",
    "",
    point -> 0 <point.a < -point.b ,
    point::SLPoint -> Inverse(point),
    CoP::CoPrime -> CoPrimesEliminatorAndPermutions(CoP),
    X -> 1,
    (x , X, N) -> x/N
    )

PlotWeylsCriterium(10^5,
    1, 
    "",
    "",
    "",
    "",
    point -> 0 <point.a < -point.b ,
    point::SLPoint -> (point),
    CoP::CoPrime -> CoPrimesEliminatorAndPermutions(CoP),
    X -> 1,
    (x , X, N) -> x/N
    )
=#

#=
PlotWeylsCriterium(10^5,
    1, 
    "",
    "",
    "",
    "",
    point -> issquare(point.a^2 + point.b^2-2) ,
    point::SLPoint -> Inverse(point),
    CoP::CoPrime -> CoPrimesEliminatorAndPermutions(CoP),
    X -> 1,
    (x , X, N) -> x/N
    )
=#

#=
PlotGpElements(10^5,"issquare(Norm(point)-2)", issquare(Norm(point)-2), point -> point)
PlotWeylsCriterium(10^6,
    1, 
    "",
    "",
    "",
    "",
    point -> issquare(Norm(point)-2) ,
    point::SLPoint -> Inverse(point),
    CoP::CoPrime -> CoPrimesEliminatorAndPermutions(CoP),
    X -> 1,
    (x , X, N) -> x/N
    )

PlotWeylsCriterium(10^5,
    1, 
    "TEST",
    "",
    "",
    "",
    point -> issquare(Norm(point)-2) ,
    point::SLPoint -> Inverse(point),
    CoP::CoPrime -> CoPrimesEliminatorAndPermutions(CoP),
    X -> 1,
    (x , X, N) -> x/N
    )
=#


#PlotGpElements(10^3.5,"SL(2,ℤ)", _ -> true, point -> point, Permutions)

#PlotGpElements(10^5,"SL(2,ℤ)∩{a²+b²-2 = n²},", point -> issquare(point.a^2 + point.b^2-2), point -> point, Permutions )
#"issquare(a^2 + b^2-2)"

#PlotGpElements(10^5,"SL(2,ℤ)∩{a²+b²-4 = n²},", point -> issquare(point.a^2 + point.b^2-4), point -> point, Permutions )

#=
PlotWeylsCriterium(10^5,
    1, 
    "Weyl's Criterium for γ ∈ SL(2,ℤ)",
    "",
    "X \n",
    L"\sum_{‖γ‖_{F^2} ≤ X} e_1(arg(γ(i)))",
    point -> true ,
    point::SLPoint -> Inverse(point),
    CoP::CoPrime -> CoPrimesEliminatorAndPermutions(CoP),
    X -> 0,
    (x , X, N) -> x/N
    )
=#
    #{} e₁(arg(γ(i)))

#=
    PlotWeylsCriterium(10^3,
    1, 
    "",
    "",
    "",
    "",
    point -> true,
    point::SLPoint -> Inverse(point),
    CoP::CoPrime -> CoPrimesEliminatorAndPermutions(CoP),
    X -> 1,
    (x , X, N) -> x/N
    )
    PlotWeylsCriterium(10^3,
    0, 
    "Weyl's Criterium for γ ∈ SL(2,ℤ)",
    "",
    "X \n",
    L"\dfrac{1}{P(X)} \sum_{‖γ‖_{F^2} ≤ X} e_1(arg(γ(i)))",
    point -> true ,
    point::SLPoint -> Inverse(point),
    CoP::CoPrime -> CoPrimesEliminatorAndPermutions(CoP),
    X -> 0,
    (x , X, N) -> x/6
    )
    =#

#=
PlotWeylsCriterium(10^4,
    0, 
    "Hyperbolic Lattice Counting problem for SL(2,ℤ)",
    "",
    "X \n",
    L"P_{SL(2,ℤ)}(X) - 6X)/X^{1/2} ",
    point -> true ,
    point::SLPoint -> Inverse(point),
    CoP::CoPrime -> Permutions(CoP),
    X -> 0,
    (x , X, N) -> (x-6*X)/(X^(1/2))
    )
=#
#=
PlotWeylsCriterium(10^8,
    0, 
    "Hyperbolic Prime for SL(2,ℤ)",
    "",
    "X \n",
    L"π_{SL(2,ℤ)}(X)(Log(X)/X) ",
    point -> isprime(Norm(point)) ,
    point::SLPoint -> Inverse(point),
    CoP::CoPrime -> Permutions(CoP),
    X -> 0,
    (x , X, N) -> (x*log(X)/X)
    )
=#
  #=  
PlotWeylsCriterium(10^4,
    0, 
    "Hyperbolic Lattice Counting problem for SL(2,ℤ)",
    "",
    "X \n",
    L"P_{SL(2,ℤ)}(X) - 6X)/X^{1/2} ",
    point -> true ,
    point::SLPoint -> Inverse(point),
    CoP::CoPrime -> Permutions(CoP),
    X -> 0,
    (x , X, N) -> (x-6*X)/(X^(1/2))
    )
    =#

    
    #=

    PlotWeylsCriterium(10^3,
    1, 
    "Weyl's Criterium for γ ∈ SL(2,ℤ) and ‖γ‖-2 is prime ",
    "",
    "X \n",
    L"\sum_{‖γ‖_{F^2} ≤ X} e_1(arg(γ(i)))",
    point -> isprime(Norm(point)-2) ,
    point::SLPoint -> point,
    CoP::CoPrime -> Permutions(CoP),
    X -> 0,
    (x , X, N) -> x/N
    )

    PlotWeylsCriterium(10^3,
    1, 
    "Weyl's Criterium for γ ∈ SL(2,ℤ) and ‖γ‖-2 is prime ",
    "",
    "X \n",
    L"\sum_{‖γ‖_{F^2} ≤ X} e_1(arg(γ(i)))",
    point -> isprime(Norm(point)-2) ,
    point::SLPoint -> point,
    CoP::CoPrime -> Permutions(CoP),
    X -> 0,
    (x , X, N) -> x/N
    )
    =#

#=
PlotWeylsCriterium(10^3,
    1, 
    "Weyl's Criterium for γ ∈ SL(2,ℤ) and ‖γ‖ is prime ",
    "",
    "X \n",
    L"\dfrac{1}{P(X)} \sum_{‖γ‖_{F^2} ≤ X} e_1(arg(γ(i)))",
    point -> isprime(Norm(point)) ,
    point::SLPoint -> point,
    CoP::CoPrime -> Permutions(CoP),
    X -> 0,
    (x , X, N) -> x/N
    ) =#
#=
    PlotWeylsCriterium(10^5,
    1, 
    "Weyl's Criterium for γ ∈ SL(2,ℤ) and ‖γ‖ is prime ",
    "",
    "X \n",
    L"\dfrac{1}{P(X)} \sum_{‖γ‖_{F^2} ≤ X} e_1(arg(γ(i)))",
    point -> issquare(point.a^2 + point.b^2+2) ,
    point::SLPoint -> point,
    CoP::CoPrime -> Permutions(CoP),
    X -> 0,
    (x , X, N) -> x/N
    )
    

    PlotWeylsCriterium(10^5,
    1, 
    "Weyl's Criterium for γ ∈ SL(2,ℤ) and ‖γ‖ is prime ",
    "",
    "X \n",
    L"\dfrac{1}{P(X)} \sum_{‖γ‖_{F^2} ≤ X} e_1(arg(γ(i)))",
    point -> isprime(Norm(point)+2) ,
    point::SLPoint -> point,
    CoP::CoPrime -> Permutions(CoP),
    X -> 0,
    (x , X, N) -> x/N
    )
    =#

    PlotGpElements(10^5,"SL(2,ℤ) and ‖γ‖+2 is a square", point -> issquare(Norm(point)+2) , point -> point)
