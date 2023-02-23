using Primes

# <Summary> 
# This Code generates files where we can count SL and more.
# Generates SL Points less and some norm.
# CountByCreateFilesNorm(10^5)
#</Summary>

# <Summary> 2x2 Matric - called SLPoint </Summary>
struct SLPoint
    a::Int
    b::Int
    c::Int
    d::Int
end

# <Summary>Give the Determinant of the SLPoint </Summary>
# <param name="point"> SLPoint</param>
# <return>int</return>
function Determinant(point::SLPoint)
    a::Int128 = point.a
    b::Int128 = point.b
    c::Int128 = point.c
    d::Int128 = point.d
    return a*d -b*c 
 end

# <Summary> Calculate the Fabenius Norm </Summary>
# <param name="point"> SLPoint </param>
# <return>int</return>
function Norm(point::SLPoint)
    a::Int128 = point.a
    b::Int128 = point.b
    c::Int128 = point.c
    d::Int128 = point.d
    return a^2 + b^2 + c^2 + d^2
end


# <Summary> Take the exp of a number  </Summary>
# <param name="base">Number - used in sums Iterators</param>
# <param name="input">Normaly the argument of a Complex number</param>
# <param name="b1">int - For all b</param>
# <param name="b2">int - For all b</param>
# <return>Int </return>
function NumberExp(base::ComplexF64, input::Vector{Float64}, b::Vector{Int}=[1])
    res = 1 + 0im
    for i in range(1,length(input))  
         res = res*exp(input[i]*1im*b[i])
    end
    base = res
end


# <Summary>Isometries return a list of the isometric Permutations of a Copirme </Summary>
# <param name="point">SLPoint</param>
# <return>A list of Isometries of SLPoint</return>

function Isometries(point::SLPoint) 
    a = point.a
    b = point.b
    c = point.c
    d = point.d

    #We still get all 8 the Isometries, so this is for that we do not get dublicates 
    if a == 1 == b
        return unique([
            SLPoint(a,b,c,d),
            SLPoint(-a,-b,-c,-d),
            SLPoint(-a,b,c,-d),
            SLPoint(a,-b,-c,d)
        ])
    end

    #We still get all 8 the Isometries, so this is for that we do not get dublicates 
    if a == 0 && b == 1
        return unique([
            SLPoint(a,b,c,d),
            SLPoint(-a,-b,-c,-d),
            SLPoint(b,a,-d,-c),
            SLPoint(-b,a,-d,c)
        ])
    end
    
    return unique([
        SLPoint(a,b,c,d),
        SLPoint(-a,-b,-c,-d),
        SLPoint(-a,b,c,-d),
        SLPoint(a,-b,-c,d),
        SLPoint(b,a,-d,-c),
        SLPoint(-b,-a,d,c),
        SLPoint(-b,a,-d,c),
        SLPoint(b,-a,d,-c)
    ])
end

# <Summary>Mobius Transformation / Action </Summary>
# <param name="point"> SLPoint </param>
# <param name="z"> Complex Number </param>
# <return>Complex Number</return>
function MobiusAction(point::SLPoint, z =1im)
    return (point.a*z + point.b)/(point.c*z + point.d)
end

# <Summary>Inverse Matrix of the SL Point</Summary>
# <param name="point"> SLPoint</param>
# <return>SLPoint</return>
function Inverse(point::SLPoint)
    return SLPoint(point.d,-point.b,-point.c,point.a)
end

# <Summary>Cayley Transformation </Summary>
# <param name="point"> SLPoint</param>
# <return>Complex Number</return>
function CayleyTransformation(z::Complex)
    return (z- 1im )/(z + 1im)
end

# <Summary>Argeumt Of Complex of the SL Point </Summary>
# <param name="point"> SLPoint</param>
# <return>int</return>
function ArgeumtOfSLPoint(point::SLPoint)
    z = CayleyTransformation(MobiusAction(point))
    iz = CayleyTransformation(MobiusAction(Inverse(point)))
    return (angle(z), angle(iz))
end

# <Summary>Quadratic Solutions k in a^2 + b^2 + (k+c)^2 + (k+d)^2 ≤ M </Summary>
# <param name="a"> a in SLPoint</param>¨
# <param name="b"> b in SLPoint</param>
# <param name="c"> c in SLPoint</param>
# <param name="d"> d in SLPoint</param>
# <param name="M"> MAX </param>
# <return>range of int</return>

function QuadraticSolutions(a,b,c,d,M)
    A = convert(Int128, a)
    B = convert(Int128, b)
    C = convert(Int128, c)
    D = convert(Int128, d)

    Mdet = (-A^4 - 2*A^2*B^2 + A^2*C^2 + 2*A*B*C*D - B^4 + B^2*D^2 + M*A^2 + M*B^2)
    if Mdet <= 0
        return 0:0
    end
    #Quadratic Solutions
    Mres1 = (A*C + B*D + sqrt(Mdet))/(A^2 + B^2)
    Mres2 = (A*C + B*D - sqrt(Mdet))/(A^2 + B^2)
    
    return floor(Int,min(Mres1, Mres2)) : ceil(Int,max(Mres1, Mres2))
end

# <Summary>Counts for a given coprime a,b and and max </Summary>
# <param name="a"> SLPoint</param>
# <param name="b"> SLPoint</param>
# <return>(int, int, ComplexF64, ComplexF64) - Count, Prime Count, WeylsSL2, WeylsC2dim )</return>
function CountBezoutGpElemtents(a,b,MAX) 
    SLTESTlist = Vector{SLPoint}()
    
    Count = 0
    PrimeCount = 0
    WeylsSL2::Number = 0.0+0.0im
    WeylsC2dim::Number = 0.0+0.0im

    gcd, d0, c0 = gcdx(a, -b)
    if gcd == 1
        c = c0
        d = d0
        for n in QuadraticSolutions(a,b,c0,d0,MAX)
            d = d0 - n * b
            c = c0 - n * a
            NewSLPoint = SLPoint(a,b,c,d)
            if Determinant(NewSLPoint) == 1 && Norm(NewSLPoint) <= MAX
                for point in Isometries(NewSLPoint)
                    Count += 1
                    if isprime(Norm(NewSLPoint))
                        PrimeCount += 1
                    end
                    arg = ArgeumtOfSLPoint(NewSLPoint)
                    WeylsSL2 = NumberExp(WeylsSL2, [arg[1]],[1])
                    WeylsC2dim = NumberExp(WeylsC2dim, [arg[1],arg[2]],[1,1])
                end
            end
        end
    end
    return (Count, PrimeCount, WeylsSL2, WeylsC2dim )
end

# <Summary>Counts SL points up to norm </Summary>
# <param name="point"> SLPoint</param>
# <return>(int, int, ComplexF64, ComplexF64) - Count, Prime Count, WeylsSL2, WeylsC2dim )</return>
function CountByCreateFilesNorm(MAX)
    m = floor(Int64,log10(MAX))
    n = 0

    Count = 0
    PrimeCount = 0
    WeylsSL2::Number = 0.0+0.0im
    WeylsC2dim::Number = 0.0+0.0im

    SLTEST2list = Vector{SLPoint}()
    res = CountBezoutGpElemtents(0,1,MAX)
    (Count, PrimeCount, WeylsSL2, WeylsC2dim ) = 
            (res[1]+Count, res[2]+PrimeCount, res[3]+WeylsSL2, res[4]+WeylsC2dim )
    res = CountBezoutGpElemtents(1,1,MAX)
    (Count, PrimeCount, WeylsSL2, WeylsC2dim ) = 
            (res[1]+Count, res[2]+PrimeCount, res[3]+WeylsSL2, res[4]+WeylsC2dim )


    println("New")

        open("CoPrimes/$m.0.dat", "w+") do file
            write(file, [1,2])
            res = CountBezoutGpElemtents(1,2,MAX)
            (Count, PrimeCount, WeylsSL2, WeylsC2dim ) = 
            (res[1]+Count, res[2]+PrimeCount, res[3]+WeylsSL2, res[4]+WeylsC2dim ) 
            write(file, [1,3]) 
            
            res = CountBezoutGpElemtents(1,3,MAX)
            (Count, PrimeCount, WeylsSL2, WeylsC2dim ) = 
            (res[1]+Count, res[2]+PrimeCount, res[3]+WeylsSL2, res[4]+WeylsC2dim )
        end    

        while filesize("CoPrimes/$m.$n.dat") > 0
            
            open("CoPrimes/$m.$n.dat", "r+") do file1
                open("CoPrimes/$m.$(n+1).dat", "w+") do file2
                    y = Array{Int64}(undef, 1, 2);
                    while ! eof(file1) 
                        read!(file1,y )
                        if 2*y[2]-y[1] <= sqrt(MAX - y[2]^2 )
                            write(file2, [y[2] , 2*y[2]-y[1]]) # 2*Co.m-Co.n,Co.m
                            res= CountBezoutGpElemtents(y[2],2*y[2]-y[1],MAX)
                            (Count, PrimeCount, WeylsSL2, WeylsC2dim ) = 
                            (res[1]+Count, res[2]+PrimeCount, res[3]+WeylsSL2, res[4]+WeylsC2dim )
                        end

                        if 2*y[2]+y[1] <= sqrt(MAX - y[2]^2 )
                            write(file2, [y[2] , 2*y[2]+y[1]]) # 2*Co.m+Co.n,Co.m
                            res = CountBezoutGpElemtents(y[2],2*y[2]+y[1],MAX)
                            (Count, PrimeCount, WeylsSL2, WeylsC2dim ) = 
                            (res[1]+Count, res[2]+PrimeCount, res[3]+WeylsSL2, res[4]+WeylsC2dim )
                        end

                        if y[2]+2*y[1] <= sqrt(MAX - y[1]^2)
                            write(file2, [y[1] , y[2]+2*y[1]]) # Co.m+2*Co.n,Co.n
                            res = CountBezoutGpElemtents(y[1],y[2]+2*y[1],MAX)
                            (Count, PrimeCount, WeylsSL2, WeylsC2dim ) = 
                            (res[1]+Count, res[2]+PrimeCount, res[3]+WeylsSL2, res[4]+WeylsC2dim )
                        end
                    end 
                end
            end
        rm("CoPrimes/$m.$n.dat")
        n += 1
    end
    rm("CoPrimes/$m.$n.dat")

    return (Count, PrimeCount, WeylsSL2, WeylsC2dim )
end

CountByCreateFilesNorm(10^9)
