#=
General model: Model parameterized for one location (Germany)
- Julia version: 1.6.3
- Authors: Kristan A. Schneider, H. Christian T. Obama, Nessma Adil M. Y., Pierre M. Ngougoue N.,
- Date created: 2021-09-06
- Date last modified: 2021-11-21
=#

using Pkg
using Plots   ### if not instaled run: import Pkg; Pkg.add("Plots")
using DifferentialEquations ### if not instaled run: import Pkg; Pkg.add("DifferentialEquations")
using LinearAlgebra   ### if not instaled run: import Pkg; Pkg.add("LinearAlgebra")
using CSV  ### if not instaled run: import Pkg; Pkg.add("CSV")
using DataFrames  ### if not instaled run: import Pkg; Pkg.add("DataFrames")


# Simulation time
tmax = 850

# Number of locations
r = 1

# Number of vaccines
V = 3

# Number of aga strata
s = 4

# Number of variants
M = 3

# Matrix of sub-population sizes r x s matrix (rows: location), (columns: age groups)
PopSize =[3969138 11365436 43730684 24722130]

# Total Population size
N = sum(PopSize) 

# Initial infections vector of length r
ininf=[[ 0. ; 10. ; 50.0 ; 5.]]   

# Erlang states
nE = [5 5 5; 5 5 5; 5 5 5; 5 5 5]  # Matrix of latent Erlang states per age group (row) and viral variant (col)
nP = [5 5 5; 5 5 5; 5 5 5; 5 5 5]  # Matrix of prodromal Erlang states per age group (row) and viral variant (col)
nI = [5 5 5; 5 5 5; 5 5 5; 5 5 5]  # Matrix of fully infectious Erlang states per age group (row) and viral variant (col)
nL = [5 5 5; 5 5 5; 5 5 5; 5 5 5]  # Matrix of late infectious Erlang states per age group (row) and viral variant (col)

# Durations of disease stages
DE = [3.5 3.5 3.5; 3.5 3.5 3.5; 3.5 3.5 3.5; 3.5 3.5 3.5]  # Duration in latent phase per age group (row) and viral variant (col)
DP = [1.0 1.0 1.0; 1.0 1.0 1.0; 1.0 1.0 1.0; 1.0 1.0 1.0]  # Duration in prodromal phase per age group (row) and viral variant (col)
DI = [5.0 5.0 5.0; 5.0 5.0 5.0; 5.0 5.0 5.0; 5.0 5.0 5.0]  # Duration in infectious phase per age group (row) and viral variant (col)
DL = [5.0 5.0 5.0; 5.0 5.0 5.0; 5.0 5.0 5.0; 5.0 5.0 5.0]  # Duration in late infectious phase per age group (row) and viral variant (col)

# Infectiousness at disease stages
cP = [0.5, 0.5, 0.5]  # Contagiousness for prodromal per viral variant
cI = [1.0, 1.0, 1.0]  # Contagiousness for fully infectious per viral variant
cL = [0.5, 0.5, 0.5]  # Contagiousness for late infectious per viral variant

# Basic reproduction number
R0     = [3., 3*1.3, 3*1.3*1.4]   # Vector of size M
tR0max = 300.0                    # Time at which the seasonal fluctuation of R0 reaches its maximum
Amp    = 0.43                     # Seasonal amplitude
lamex  = 30.0                     # External force of infection (this will be split to the age groups in the respective locations according to their contact behavior)

# Time eacg viral variant was introduced
mutint = [-20 290 475]

# Severity of infection
## Fraction of sick individuals per age group (row) and viral variant (col) sxM Matrix
fsick = [0.15 0.15 0.15; 0.30 0.30 0.30; 0.65 0.65 0.65; 0.70 0.70 0.70] 

## Fraction of partially immunized sick individuals per age group (row) and viral variant (col) sxM Matrix
fsickPI = [[0.10 0.11 0.10; 0.11 0.12 0.11; 0.12 0.13 0.12 ],
           [0.25 0.22 0.25; 0.26 0.23 0.26; 0.27 0.24 0.27 ],
           [0.55 0.60 0.55; 0.56 0.61 0.56; 0.57 0.62 0.57 ],
           [0.60 0.60 0.60; 0.60 0.60 0.60; 0.60 0.60 0.60 ]] 

## Fraction of dead individuals per age group (row) and viral variant (col) sxM Matrix
fdead = [0.0001 0.0002 0.0004; 0.0001 0.0002 0.0004; 0.0010 0.0011 0.0012; 0.24 0.24 0.24]

## Fraction of dead partially immunized individuals per age group (vector-element): entry is MxV Matrices per viral variant (row) and vaccine (col))
fdeadPI = [[0.0        0.0      0.0     ; 0.0       0.0      0.0     ;0.0        0.0      0.00],
           [5.0e-6     5.5e-6   5.0e-6  ; 1.0e-5    1.6e-5   1.0e-5  ;2.0e-5     2.5e-5   2.0e-5],
           [0.000055   0.00006  0.000055; 0.000055  0.00006  0.000055;0.000055   0.00006  0.000055],
           [0.004      0.0045   0.004   ; 0.004     0.0045   0.004   ;0.004      0.0045   0.004 ]]

# Fraction of sick individuals that die
fsickdead = fsick .* fdead

# Fraction of partially immunized sick individuals that die per age group
fsickdeadPI = fsickPI
for a = 1:s
    fsickdeadPI[a] = fsickPI[a] .* fdeadPI[a]
end

# Vaccination onset
# Array of dimension V (number of vaccines prducts) x s (age strata) x l (locations) describes when the VC start, per age group (row) and location (col)
tvax= [[ 850  490  400  310 ],
        [850  850  400  360 ],
        [850  850  430  430 ]]

## Waiting time for vaccination
## Vector of size V (number of vaccines) - in each entry r x s matrix with location (row) and age group (col)
## If vaccination rate changes for any vaccine, location, or age group the vector needs to be redefined
Dvax1= [[ 1  1  1  280 ],
        [ 1  1  1  1 ],
        [ 1  1  1  1 ]]

Dvax2= [[ 1  1  1  160 ],
        [ 1  1  1  400 ],
        [ 1  1  1  1 ]]

Dvax3= [[ 1 1 240 160 ],
        [ 1 1 240 400 ],
        [ 1 1 1   1 ]]

Dvax4= [[ 1  1  240  140 ],
        [ 1  1  240  340 ],
        [ 1  1  240  400 ]]

Dvax5= [[ 1  170  240  140 ],
        [ 1  1    240  340 ],
        [ 1  1    240  400 ]]

## Vector of waiting times for vaccination
Dvaxlist = [Dvax1, Dvax2, Dvax3, Dvax4, Dvax5, Dvax5] 

## Times vaccination rates change corresponding to matrices Dvax1, ...,Dvax5
tvaxchange = [360, 400, 430, 490, 590]

# Vector with proportion of vaccinable individuals per location and entries being vectors of length s (age groups)
propvax = [[0.0 ; 0.36 ; 0.85 ; 0.95 ]]   

## For mandatory vaccine activate if condition
if 1==0
    ## Waiting time for vaccination
    ## Vector of size V (number of vaccines) - in each entry r x s matrix with location (row) and age group (col)
    ## If vaccination rate changes for any vaccine, location, or age group the vector needs to be redefined

    Dvax1= [[ 1  1  1  300 ],
            [ 1  1  1  1 ],
            [ 1  1  1  1 ]]

    Dvax2= [[ 1  1  1  180 ],
            [ 1  1  1  600 ],
            [ 1  1  1  1 ]]

    Dvax3= [[ 1 1 280 160 ],
            [ 1 1 300 450 ],
            [ 1 1 1   1 ]]

    Dvax4= [[ 1  1  240  120 ],
            [ 1  1  300  340 ],
            [ 1  1  300  400 ]]

    Dvax5= [[ 1  300  270  140 ],
            [ 1  1    300  340 ],
            [ 1  1    300  400 ]]

    Dvax6= [[ 1  30  30  30 ],
                [ 1  1    240  340 ],
                [ 1  1    200  100 ]]

    ## Vector of waiting times for vaccination
    Dvaxlist = [Dvax1, Dvax2, Dvax3, Dvax4, Dvax5, Dvax6, Dvax6]

    ## Times vaccination rates change corresponding to matrices Dvax1, ...,Dvax6
    tvaxchange = [360, 400, 430, 490,  636, 700]

    # Vector with proportion of vaccinable individuals per location and entries being vectors of length s (age groups)
    propvax = [[0.0 ; 0.51 ; 0.95 ; 0.97 ]]   
end

# Vaccine effect for each vaccine: fS [fS(Im), fS(PI), fS(NI)], fE [fE(Inf), fE(PI)], fE(NI)] , fP, ... as vector of length V with entries being 5 x 3 matrices (rows stages, columns IM, PI, NI)
f= [[0.25 0.7 0.05;  0. 0.5 0.5;  0. 0.5 0.5;  0. 0.5 0.5;  0. 0.5 0.5],
    [0.3 0.65 0.05;  0. 0.5 0.5;  0. 0.5 0.5;  0. 0.5 0.5;  0. 0.5 0.5],
    [0.3 0.65 0.05;  0. 0.5 0.5;  0. 0.5 0.5;  0. 0.5 0.5;  0. 0.5 0.5]]

# Waiting time for vaccine to immunize for age strata (rows) and vaccines (cols)
DA = [30 50 15; 30 50 15; 30 50 15; 30 50 15]

# Vaccine protection from variants (MxV Matrix with viral variants (row), vaccine (col)
g = [0.5 0.5 0.55; 0.60 0.60 0.60; 0.60 0.6 0.6] 
h = [0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0]

# Effect of partial immunity on transmission
## During prodromal phase m x v matrix
pPmv = [ 0.66 0.70 0.66
         0.66 0.70 0.66
         0.66 0.70 0.66]

## During fully infectious phase m X v matrix
pImv = [ 0.66 0.70 0.66
         0.66 0.70 0.66
         0.66 0.70 0.66]

## During late infectious phase m X v matrix
pLmv = [ 0.66 0.70 0.66
         0.66 0.70 0.66
         0.66 0.70 0.66]

#############################
# Case isolation
#############################

Qmax  = [200]    # Vector of length r: maximum capacity per location per 100.000
fiso  = [0.58]   # Probability of being isolated per location, vector of length r
tiso1 = [10.]    # Time isolation starts
tiso2 = [tmax]   # Time isolation ends
phome = 0.75     # Contacts prevented at home isolation

#############################
# General contact reduction
## The contact matrix and the contact reductions need to be implemented to generate the contact matrices X(t),
## this can be done individuzally and the code hence needs to be adjusted properly.
##t0 = 0.0   -> Feb. 25, 2020 - First reported case in Germany.
#############################

## Times for weather adjustment
tw1 = 380.0    # Mar. 11, 2021
tw2 = 430.0    # Apr. 30, 2021
tw3 = 520.0    # Jul. 29, 2021
tw4 = 540.0    # Aug. 08, 2021

## Weather adjustment
w_adj = [20, -30]

## Times contact reductions change (NN such time points)
t1  = 40.0    # Apr. 02, 2020 - First lockdown
t2  = 85.0    # May  20, 2020 - Relaxation
t3  = 97.0    # Jun. 01, 2020 - summer vacations
t4  = 170.0   # Aug. 23, 2020 - end summer vacations
t5  = 190.0   # Sep. 02, 2020 - week measures
t6  = 245.0   # Oct. 27, 2020 - soft lockdown
t7  = 280.0   # Dec. 01, 2020 - hard christmas lockdown
t8  = 303.0   # Dec. 24, 2020 -
t9  = 355.0   # Feb. 14, 2021 -
t10 = 425.0   # Apr. 25, 2021 - emergency break
t11 = 490.0   # Jun. 29, 2021 - End emergency break
t12 = 540.0   # Aug. 23, 2021 - 3G rule (Sachsen)
t13 = 621.0   # Nov. 07, 2021 - 2G rule
t14 = 636.0   # Nov. 22, 2021 - hypothetical shcool closures
t15 = tmax    # End of simulation

## Interval of school vacation (summer vacation)
tschool = [501.0 537.0; 866.0 904.0; 1100.0 1150.0]

## Time Emergency Break (EB) starts
tEB = t10

## Time emergency break ends
tEBstop = t11

## Time Emergency Break (EB) starts again (use for hypothetical roll back of EB)
tEB2 = tmax

## Time emergency break stops again
tEBstop2 = tmax

## Breaking time of general contact reductions
conttime = copy([0, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15])

## Contact reduction per break point for Home, Others, School, Work
## NN x 6 matrix with amount of contacts being reduced Home, Home_old (oldest age group), Others, Others_old (Oldest age group), School, Work at various time points

ContRed = [ 0.35  0.05  0.95  0.20  0.95  0.70;    #t1:  Apr. 03, 2020 - First lockdown
            0.15  0.05  0.60  0.20  0.75  0.50;    #t2:  May  20, 2020 - Relexation
            0.00  0.00  0.40  0.10  1.00  0.20;    #t3:  Jun. 01, 2020 - summer vacations
            0.10  0.00  0.50  0.10  0.50  0.40;    #t4:  Aug. 13, 2020 - end summer vctations
            0.15  0.05  0.70  0.20  0.75  0.50;    #t5:  Sep. 02, 2020 - week measures
            0.33  0.30  0.80  0.70  1.00  0.75;    #t6:  Oct. 27, 2020 - soft lockdown
            0.35  0.30  0.80  0.75  1.00  0.75;    #t7:  Dec. 01, 2020 - hard christmas lockdown
            0.30  0.35  0.95  0.75  1.00  0.85;    #t8:  Dec. 24, 2020
            0.25  0.20  0.75  0.60  0.75  0.60;    #t9:  Feb. 14, 2021
            0.00  0.00  0.00  0.00  0.00  0.00;    #t10: Apr. 25, 2021 - emergency break
            0.10  0.05  0.40  0.25  0.50  0.20;    #t11: Jun. 29, 2021 - End emergency break
            0.30  0.15  0.85  0.40  0.75  0.65;    #t12: Aug. 23, 2021 - 3G rule (Schleswig Holstein)
            0.30  0.15  0.85  0.40  0.75  0.80;    #t13: Nov. 07, 2021 - 2G rule
            0.30  0.15  0.85  0.40  0.75  0.80;    #t14: Nov. 22, 2021 - hypothetical schools closed
            0.30  0.15  0.85  0.40  0.75  0.80]    #t15: End of simulation

## Incidence threshold values (per 100,000) triggering EB contact reductions (NN1 time points)
Incid_Trig = [10.0, 50.0, 100.0, 180.0]

## Incidence based contact reductions for emergency break ((NN1+1) x 6 matrix)
ContRedInd = [ 0.00  0.00  0.50  0.10  0.25  0.30;
               0.25  0.05  0.60  0.20  0.50  0.55;
               0.30  0.10  0.85  0.60  0.75  0.65;
               0.40  0.20  0.95  0.75  0.75  0.85;
               0.40  0.35  0.95  0.75  1.00  0.85]

# Contact behavior (case of Germany)
## Contacts at home
XHome = [[ 1751573.48890175  1252883.17470613  6905537.24799685      215245.301789968;
           1252883.17470613  10707295.9471861  18242448.1345504      752354.096153554;
           6905537.24799685  18242448.1345504  51714202.4788967      3.454940099080605e6;
           215245.301789968  752354.096153554  3.454940099080605e6   0.0]]

## Contacts at home for elderlies (60+ years)
XHomeOld = [[ 0.0  0.0  0.0                  0.0; 
              0.0  0.0  0.0                  0.0;
              0.0  0.0  0.0                  3.454940099080605e6;
              0.0  0.0  3.454940099080605e6  16545972.3421171]]

# Contact at other places
XOthers = [[ 2052141.67004117  1241705.04837431  5003947.84744736    1331986.38972126;
             1241705.04837431  15595046.381503   12421035.394769     2453842.8726374;
             5003947.84744736  12421035.394769   80096635.5761702    1.04824470460359e7;
             1331986.38972126  2453842.8726374   1.04824470460359e7  0.0]]

## Contacts at other places for elderlies (60+ years)
XOthersOld = [[ 0.0  0.0  0.0                 0.0; 
                0.0  0.0  0.0                 0.0;
                0.0  0.0  0.0                 1.04824470460359e7;
                0.0  0.0  1.04824470460359e7  17359981.2818305]]

## Contacts in schools
XSchool = [[  5268910.91638032      1063592.70685887  1235765.71316362  3.08138382710115e-60;
              1063592.70685887      15472569.0475962  9573684.72254725  341276.930372909;
              1235765.71316362      9573684.72254725  6426785.46477011  57093.347800278;
              3.08138382710115e-60  341276.930372909  57093.347800278   116268.749823077]]

## cContacts at work.
XWork = [[  0.0                   0.0                0.0                5.18649431290478e-36;
            0.0                   2610993.13268279   5739939.11461415   20284.8053163813;
            0.0                   5739939.11461415   141403258.818845   709188.636019723;
            5.18649431290478e-36  20284.8053163813   709188.636019723   7125.98797403161]]

## Fraction of general contacts reduced at home, others, school, work at each breaking point
redu = [vcat([0 0 0 0 0 0],ContRed,[0 0 0 0 0 0])]

## Fraction of contacts reduced during emergency break at home, others, school, work at each breaking point
reduinz=[vcat(ContRedInd)]

# General contact reduction
## Proportion of the population moving to work
mobilwork = [1.0] 

## Proportion of the population moving to other places (vacation, etc.)
mobilother = [1.0] 

Xmat = Array{Union{Missing, Any}}(missing, length(conttime))
XmatEB = Array{Union{Missing, Any}}(missing, (length(Incid_Trig)+1)^r)
XmatNoSch= Array{Union{Missing, Any}}(missing, (length(Incid_Trig)+1)^r)

for t = 1:length(conttime)
      Contmat = zeros(s*r,s*r)
      for l in 1:r
          A = XHome[l] .* (1-redu[l][t,1])+ XHomeOld[l] .* (1-redu[l][t,2])+  XOthers[l] .* (1-redu[l][t,3]) .* mobilother[l,l]+ XOthersOld[l] .* (1-redu[l][t,4]) .* mobilother[l,l]+ XSchool[l] .* (1-redu[l][t,5])  + XWork[l] .* (1-redu[l][t,6]) .* mobilwork[l,l]
          Contmat[(l-1)*s+1:l*s,(l-1)*s+1:l*s] = N^2 ./(repeat(PopSize[l,1:s],inner=(1,s)) .* repeat(PopSize[l,1:s],inner=(1,s))') .*A
          for m in (l+1):r
              B =  XOthers[l] .* (1-redu[l][t,3]) .* mobilother[l,m] + XOthersOld[l] .* (1-redu[l][t,4]) .* mobilother[l,m] + XWork[l] .* (1-redu[l][t,6]) .* mobilwork[l,m]
              C =  XOthers[m] .* (1-redu[m][t,3]) .* mobilother[m,l] + XOthersOld[m] .* (1-redu[m][t,4]) .* mobilother[m,l] + XWork[m] .* (1-redu[m][t,6]) .* mobilwork[m,l]
              C1 = (B+C)/2
              C1 = N^2 .* C1 ./ (repeat(PopSize[l,1:s],inner=(1,s)) .* repeat(PopSize[m,1:s],inner=(1,s))')
              Contmat[(l-1)*s+1:l*s,(m-1)*s+1:m*s] = C1
              Contmat[(m-1)*s+1:m*s,(l-1)*s+1:l*s] = C1'
          end
      end
      Xmat[t] = Contmat
end

inctr = (length(Incid_Trig)+1)
AA = zeros(Int64, inctr^r,r)
AA[1:inctr,1] = 1:inctr

for l = 2:r
    AA[1:inctr^l,1:(l-1)] = repeat(AA[1:inctr^(l-1),1:(l-1)],outer=(inctr,1))
    AA[1:inctr^l,l] = repeat(1:inctr,inner=inctr^(l-1))
end

base = inctr .^(0:(r-1))

# Contact reduction for emergency break (EB)
for t = 1:inctr^r
    ind = AA[t,:]

    Contmat = zeros(s*r,s*r)
    for l in 1:r
        A = XHome[l] .* (1-reduinz[l][ind[l],1])+ XHomeOld[l] .* (1-reduinz[l][ind[l],2])+ XOthers[l] .* (1-reduinz[l][ind[l],3]) .* mobilother[l,l]+ XOthersOld[l] .* (1-reduinz[l][ind[l],4]) .* mobilother[l,l]+ XSchool[l] .* (1-reduinz[l][ind[l],5])  + XWork[l] .* (1-reduinz[l][ind[l],6]) .* mobilwork[l,l]
        Contmat[(l-1)*s+1:l*s,(l-1)*s+1:l*s] = N^2 ./(repeat(PopSize[l,1:s],inner=(1,s)) .* repeat(PopSize[l,1:s],inner=(1,s))') .*A
        for m in (l+1):r
            a = max(ind[l],ind[m])
            if a ==inctr  ### if maximum incidence is surpassed 20 km restriction area, only work related travel
                B = XWork[l] .* (1-reduinz[l][a,6]) .* mobilwork[l,m]
                C = XWork[m] .* (1-reduinz[m][a,6]) .* mobilwork[m,l]
            else
                B =  XOthers[l] .* (1-reduinz[l][a,3]) .* mobilother[l,m] + XOthersOld[l] .* (1-reduinz[l][a,4]) .* mobilother[l,m] + XWork[l] .* (1-reduinz[l][a,6]) .* mobilwork[l,m]
                C =  XOthers[m] .* (1-reduinz[m][a,3]) .* mobilother[m,l] + XOthersOld[m] .* (1-reduinz[m][a,4]) .* mobilother[m,l] + XWork[m] .* (1-reduinz[m][a,6]) .* mobilwork[m,l]
            end
            C1 = (B+C)/2
            C1 = N^2 .* C1 ./ (repeat(PopSize[l,1:s],inner=(1,s)) .* repeat(PopSize[m,1:s],inner=(1,s))')
            Contmat[(l-1)*s+1:l*s,(m-1)*s+1:m*s] = C1
            Contmat[(m-1)*s+1:m*s,(l-1)*s+1:l*s] = C1'
        end
    end
    XmatEB[t] = Contmat
end

# Contact reduction during EB assuming no school (Summer vacations)
for t = 1:inctr^r   ## no school
    ind = AA[t,:]

    Contmat = zeros(s*r,s*r)
    for l in 1:r
        A = XHome[l] .* (1-reduinz[l][ind[l],1])+ XHomeOld[l] .* (1-reduinz[l][ind[l],2])+ XOthers[l] .* (1-reduinz[l][ind[l],3]) .* mobilother[l,l]+ XOthersOld[l] .* (1-reduinz[l][ind[l],4]) .* mobilother[l,l] + XWork[l] .* (1-reduinz[l][ind[l],6]) .* mobilwork[l,l]
        Contmat[(l-1)*s+1:l*s,(l-1)*s+1:l*s] = N^2 ./(repeat(PopSize[l,1:s],inner=(1,s)) .* repeat(PopSize[l,1:s],inner=(1,s))') .*A
        for m in (l+1):r
            a = max(ind[l],ind[m])
            if a ==inctr  ### if maximum incidence is surpassed 20 km restriction area, only work related travel
                B = XWork[l] .* (1-reduinz[l][a,6]) .* mobilwork[l,m]
                C = XWork[m] .* (1-reduinz[m][a,6]) .* mobilwork[m,l]
            else
                B =  XOthers[l] .* (1-reduinz[l][a,3]) .* mobilother[l,m] + XOthersOld[l] .* (1-reduinz[l][a,4]) .* mobilother[l,m] + XWork[l] .* (1-reduinz[l][a,6]) .* mobilwork[l,m]
                C =  XOthers[m] .* (1-reduinz[m][a,3]) .* mobilother[m,l] + XOthersOld[m] .* (1-reduinz[m][a,4]) .* mobilother[m,l] + XWork[m] .* (1-reduinz[m][a,6]) .* mobilwork[m,l]
            end
            C1 = (B+C)/2
            C1 = N^2 .* C1 ./ (repeat(PopSize[l,1:s],inner=(1,s)) .* repeat(PopSize[m,1:s],inner=(1,s))')
            Contmat[(l-1)*s+1:l*s,(m-1)*s+1:m*s] = C1
            Contmat[(m-1)*s+1:m*s,(l-1)*s+1:l*s] = C1'
        end
    end
    XmatNoSch[t] = Contmat
end
#end

###############################################################
###### Derived parameters
###############################################################
# maximum capacity per 100,000 per location
Qmax = Qmax .* sum(PopSize,dims=2) ./100000  

X = copy(Xmat[1])

lamexvec = zeros(r,s)
for l = 1:r
    XX = copy(Xmat[1])[((1:s) .+(l-1)*s),((1:s) .+(l-1)*s)]
    lamexvec[l,:] = (zeros(s) .+1)*lamex[l]
end

mutinttimes = sort(unique(collect(Iterators.flatten(mutint))))

lamexmat = fill(zeros(r,s,M), length(mutinttimes))
for t=1:length(mutinttimes)
    pick= sum(mutint .<= mutinttimes[t],dims=2)

    lamexmat0=zeros(r,s,M)
    for l=1:r
            lamexmat0[l,:,pick[l]] = collect(Iterators.flatten(lamexvec[l,:] .* N ./PopSize[l,:]))
    end
    lamexmat[t]=    lamexmat0
end

lamexmat[2] =lamexmat[2] ./2
lamexmat[3] =lamexmat[3] 

# Rates Beta
cPmat = repeat(cP,inner=(1,s))'
cImat = repeat(cI,inner=(1,s))'
cLmat = repeat(cL,inner=(1,s))'

betaP = cPmat
betaI = cImat
betaL = cLmat

nEP   = nE + nP
nEPI  = nEP + nI
nEPIL = nEPI + nL

# Number of infected, dead, and recoverd compartments per age
CnumInf = sum(nEPIL .*(2*(V+1)).+2,dims=2)  

# Total Number of compartments per age
Cnum = sum(nEPIL .*(2*(V+1)).+2,dims=2) .+(3*V+2) 

# Erl[s+1] number of compartments for each location
Erl = accumulate(+,[0 transpose(Cnum)]) 

# Matrix with indices of S(U) per age group (row) per loction (col)
IndSU = repeat(transpose(collect(0:(r-1))).*Erl[s+1],outer=(s,1))
IndSU = IndSU .+ repeat(Erl[1:s].+1,outer=(1,r)) 

# Index where the first Erlang states start
Erl1 = accumulate(+,[0 transpose(Cnum)])[1:s] .+(3*V+3)
IndE1Uam = zeros(Int64,s,M,r)

for l = 1:r
    IndE1Uam[1:s,1:M,l] = hcat(zeros(Int64, s,1), accumulate(+,nEPIL .*(2*(V+1)).+2,dims=2))[:,1:M] .+ Erl1 .+ (l-1)*Erl[s+1]
end

# Arrays of transition rates for (U) compartments
rates = Array{Union{Missing, Any}}(missing, s, M)
for a = 1:s
    for m = 1:M
        eps1 = nE[a,m] .*repeat([1/DE[a,m]],nE[a,m])
        phi = nP[a,m] .*repeat([1/DP[a,m]],nP[a,m])
        gamma = nI[a,m] .*repeat([1/DI[a,m]],nI[a,m])
        delta =  nL[a,m] .*repeat([1/DL[a,m]],nL[a,m])
        rates[a,m] = vcat(eps1, phi, gamma, delta)
    end
end

# Next generation matrix
cPIL = Array{Union{Missing, Any}}(missing, s, M)
for a = 1:s
    for m = 1:M
        cEnew = repeat([0],nE[a,m])
        cPnew = repeat([cP[m]],nP[a,m])
        cInew = repeat([cI[m]],nI[a,m])
        cLnew =  repeat([cL[m]],nL[a,m])
        cPIL[a,m] = vcat(cEnew, cPnew, cInew, cLnew)
    end
end

R0new = zeros(M)

for m = 1:M
    # Total number of states for viral variant m
    nstat = nE[:,m]+nP[:,m]+nI[:,m]+nL[:,m]  
    nstat2 = accumulate(+, nstat)
    nstat1 = accumulate(+, vcat(0,nstat)).+1
    numstates = sum(nstat)
    dimV = numstates*r

    # Matrix V (transition rates)
    dV = zeros((dimV,dimV)) 

    # Matrix F (new infections)
    dF = zeros((dimV,dimV)) 
    for l = 1:r
        for a = 1:s
            ind1 = (l-1)*numstates + nstat1[a]
            ind2 = (l-1)*numstates + nstat2[a]
            ra = length(rates[a,m])
            dV[ind1:ind2,ind1:ind2] =  fill(0.,(ra,ra))+ Bidiagonal(-rates[a,m],rates[a,m][1:(ra-1)],:L)
            for lt = 1:r
                for at = 1:s
                    ind1t = (lt-1)*numstates + nstat1[at]
                    ind2t = (lt-1)*numstates + nstat2[at]
                    dF[ind1,ind1t:ind2t] = X[(l-1)*s+a,(lt-1)*s+at] .* cPIL[at,m]* PopSize[l,a] ./N
                end
            end
        end
    end
    R0new[m] = R0[m]/maximum(map.(abs,eigvals(-dF*inv(dV))))*(1+Amp*cos(-2*pi*tR0max/365))
end

# Rates for vaccine effects to manifest
alpha = 1 ./DA

# Look up table for vaccination rates
vaxtime1 =  unique(sort(collect(Iterators.flatten(tvax))))
vaxtime =unique(sort(collect(vcat(vaxtime1,tvaxchange))))

T = length(vaxtime)
vaxrates = Array{Union{Missing, Any}}(missing, T+1)
for t = 1:(T)
    tmp1 = map(x -> x .< vaxtime[t], tvax)
    idx_Dvax = sum(tvaxchange .< vaxtime[t])+1
    Dvax = Dvaxlist[idx_Dvax]
    tmp2 = map( y -> 1 ./y, Dvax)
    tmp3 = map( (x,y) -> x.*y, tmp1,tmp2)
    tmp4 =  zeros((s,r,V))
    for v= 1:V
        tmp4[:,:,v] = tmp3[v]'
    end
    vaxrates[t] = tmp4
end

tmp4 =  zeros((s,r,V))
tmp3 = map( x -> 1 ./x, Dvaxlist[length(Dvaxlist)])

for v= 1:V
    tmp4[:,:,v] = tmp3[v]'
end
vaxrates[T+1] = tmp4

fratesS = Array{Float32}(undef, (V,3)) # Here, 3 because Im, PI, NI
frates = Array{Union{Missing, Any}}(missing, V)

IndE1 = Array{Union{Missing, Any}}(missing, r)
IndLnL = Array{Union{Missing, Any}}(missing, r)
IndRInf = Array{Int32}(undef, (s, M, r))
INDvec = Array{Union{Missing, Any}}(missing, r)

# Vector of indices for the prodromal, fully infectious, and late infectious indivisuals, respectively
INDvecP = Array{Union{Missing, Any}}(missing, r) 
INDvecI = Array{Union{Missing, Any}}(missing, r)
INDvecL = Array{Union{Missing, Any}}(missing, r)

for l = 1:r
    indage = Array{Union{Missing, Any}}(missing, s)
    indP = Array{Union{Missing, Any}}(missing, s)    # Indices for Prodromals per location r, all age groups a, and all strains M.
    indI = Array{Union{Missing, Any}}(missing, s)    # Indices for I
    indL = Array{Union{Missing, Any}}(missing, s)    # Indices for L
    indE1l = Array{Union{Missing, Any}}(missing, s)  # Indices for E1
    indLnLl = Array{Union{Missing, Any}}(missing, s) # Indices for LnL

    for a = 1:s
        indm = Array{Union{Missing, Any}}(missing, M)
        indPm = Array{Union{Missing, Any}}(missing, M)
        indIm = Array{Union{Missing, Any}}(missing, M)
        indLm = Array{Union{Missing, Any}}(missing, M)
        indE1la = Array{Union{Missing, Any}}(missing, M)
        indLnLla = Array{Union{Missing, Any}}(missing, M)

        for m = 1:M
            indcomp = Array{Union{Missing, Any}}(missing, 2*V+2) # Number of compartments at a given disease state
            indcompP = Array{Union{Missing, Any}}(missing, 2*V+2)
            indcompI = Array{Union{Missing, Any}}(missing, 2*V+2)
            indcompL = Array{Union{Missing, Any}}(missing, 2*V+2)

            pick0 = IndE1Uam[a,m,l]:(IndE1Uam[a,m,l] + nEPIL[a,m]-1)    # Indices of all infected states
            pick = (IndE1Uam[a,m,l]+1):(IndE1Uam[a,m,l] + nEPIL[a,m]-1) # Indices of infected states starting from E2
            pick1 = (IndE1Uam[a,m,l]):(IndE1Uam[a,m,l] + nEPIL[a,m]-2)  # Indicesof infected states until L_nL

            indcompabc = Array{Union{Missing, Any}}(missing, 5)
            indcompabc[1] = pick0 .+ 0 .* (nEPIL[a,m])
            indcompabc[2] = pick .+ 0 .* (nEPIL[a,m])
            indcompabc[3] = pick1 .+ 0 .* (nEPIL[a,m])
            indcompabc[4] = 1:(nEPIL[a,m]-1) .*1+.+ 0 .* (nEPIL[a,m]-1)
            indcompabc[5] = 2:nEPIL[a,m] .+ 0 .* (nEPIL[a,m]-1)
            indcomp[1] = indcompabc
            indE1lam = Array{Union{Missing, Any}}(missing, 2*V+2)
            indLnLlam = Array{Union{Missing, Any}}(missing, V+2)
            indLnLlamPI = Array{Union{Missing, Any}}(missing, V)
            indE1lam[1] = IndE1Uam[a,m,l]
            indLnLlam[1] = IndE1Uam[a,m,l] + nEPIL[a,m]-1

            pickP =  (IndE1Uam[a,m,l] + nE[a,m]):(IndE1Uam[a,m,l] + nEP[a,m]-1)     # Indices of the unvaccinated prodromals PU.
            pickI =  (IndE1Uam[a,m,l] + nEP[a,m]):(IndE1Uam[a,m,l] + nEPI[a,m]-1)   # Indices of the unvaccinated fully infectious IU.
            pickL =  (IndE1Uam[a,m,l] + nEPI[a,m]):(IndE1Uam[a,m,l] + nEPIL[a,m]-1) # Indices of the the unvaccinated late infectious LU.

            indcompP[1] =  pickP # Indices of the unvaccinated prodromals PU.
            indcompI[1] =  pickI # Indices of the unvaccinated fully infectious IU.
            indcompL[1] =  pickL # Indices of the the unvaccinated late infectious LU.

            for v in 2:(2*V+2)
                indcompabc = Array{Union{Missing, Any}}(missing, 3)
                indcompabc[1] = pick0 .+ (v-1) .* (nEPIL[a,m])
                indcompabc[2] = pick .+ (v-1) .* (nEPIL[a,m])
                indcompabc[3] = pick1 .+ (v-1) .* (nEPIL[a,m])
                indcomp[v] = indcompabc
                indE1lam[v] = IndE1Uam[a,m,l] .+ (v-1) .* (nEPIL[a,m])
                indcompP[v] =  pickP .+ (v-1) .* (nEPIL[a,m])  # Indices of the vaccinated prodromals PV.
                indcompI[v] =  pickI .+ (v-1) .* (nEPIL[a,m])  # Indices of the vaccinated fully infectious IV.
                indcompL[v] =  pickL .+ (v-1) .* (nEPIL[a,m])  # Indices of the vaccinated late infectious LV per vaccine
            end

            for v in 2:(V+1)
                indLnLlam[v] = indE1lam[v] + nEPIL[a,m] -1
            end

            indLnLlam[V+2] =   IndE1Uam[a,m,l] .+ (2*V+2) .* (nEPIL[a,m]) - 1

            for v in 1:V 
                indLnLlamPI[v] = indE1lam[v+V+1] + nEPIL[a,m] -1
            end

            IndRInf[a,m,l] =  IndE1Uam[a,m,l] + (2*V+2) * nEPIL[a,m]
            indm[m] = indcomp
            indE1la[m] = indE1lam
            indLnLla[m] = [indLnLlam, indLnLlamPI]

            indPm[m] = indcompP
            indIm[m] = indcompI
            indLm[m] = indcompL
        end
        indage[a] = indm      # Indices per age group
        indE1l[a] = indE1la
        indLnLl[a] = indLnLla

        indP[a] = indPm       # Indices of prodromals per age
        indI[a] = indIm       # Indices of fully infectious per age
        indL[a] = indLm       # Indices of late infectious per age
    end
    INDvec[l] = indage
    IndE1[l] = indE1l
    IndLnL[l] = indLnLl

    INDvecP[l] = indP  
    INDvecI[l] = indI  
    INDvecL[l] = indL
end

for v = 1:V
    fratesS[v,:] = f[v][1,:]
    rates1 = Array{Union{Missing, Any}}(missing, s, M)
    for a = 1:s
        for m = 1:M
            # Transtition rates for latent, prodromal, fully infectious, and late infectious individuals respectively.
            afE = repeat(transpose(f[v][2,:] ./DA[a,v]),outer=(nE[a,m],1)) 
            afP = repeat(transpose(f[v][3,:] ./DA[a,v]),outer=(nP[a,m],1))
            afI = repeat(transpose(f[v][4,:] ./DA[a,v]),outer=(nI[a,m],1))
            afL = repeat(transpose(f[v][5,:] ./DA[a,v]),outer=(nL[a,m],1))
            rates1[a,m] = vcat(afE, afP, afI, afL)
        end
    end
    frates[v] = rates1
end

## build matrix with values f for E1 compartments
fE1 =  Array{Float32}(undef, (M,V)) # f3 =  Im, PI, NI
for m = 1:M
    for v = 1:V
        fE1[m,v] = f[v][2,2]
    end
end

# Death rate
delta = nL ./DL

## Store model parameters
IND = [IndSU, IndE1Uam, IndE1, nEPIL, INDvec, IndLnL, IndRInf] # Indices for unvaccinated Susceptibles, Latents...
INDInf = [INDvecP, INDvecI, INDvecL]                           # Indices for infective individuals
RATES = [rates, alpha, vaxrates, vaxtime, fratesS,  frates, fE1, fsickdead, delta, betaP, betaI, betaL] # rates
tspan=(0.,tmax)
thres = 1:1.:tmax                                              # Callback incidence time points
pars = [N, g, h, R0new,  Amp, tR0max]
plam = [INDInf]                                                # Parameters for the force of infection lambda
Idx_Incd = IndRInf[s,M,r]+1
p = [IND, RATES, pars, [Xmat, XmatEB, XmatNoSch, conttime, tEB, tEBstop, tEB2, tEBstop2], INDInf, lamexmat, thres, Idx_Incd, r, X, N, Incid_Trig, tschool, mutint,[pPmv,pImv,pLmv]]

##########################################################################
###### Main functions force of infection + implementation of ODE system
##########################################################################

## The force of infection
function lambda(u,INDInf,Qmax,fiso,betaP,betaI,betaL,R0new,Amp,X,tR0max,lamexmat,t,mutint,pPmv,pImv,pLmv, mutinttimes)

    EffInf=Array{Float32}(undef, (r*s,M)) # Effective number of infectious individuals per location age group and viral variant

    PsumU=Array{Float32}(undef, (r,s,M))  # Sum of the infective prodromal individuals
    IsumU=Array{Float32}(undef, (r,s,M))  # Sum of the infective fully infectious individuals
    LsumU=Array{Float32}(undef, (r,s,M))  # Sum of the infective late infectious individuals

    PsumPIv=Array{Float32}(undef, (r,s,M,V)) # Sum of the infective prodromal individuals, Partially imunized inds stratified by vaccine
    IsumPIv=Array{Float32}(undef, (r,s,M,V)) # Sum of the infective fully infectious individuals, Partially imunized inds stratified by vaccine
    LsumPIv=Array{Float32}(undef, (r,s,M,V)) # Sum of the infective fully infectious individuals, Partially imunized inds stratified by vaccine

    PsickU=Array{Float32}(undef, (r,s,M)) # Sum of the infective prodromal individuals
    IsickU=Array{Float32}(undef, (r,s,M)) # Sum of the infective fully infectious individuals
    LsickU=Array{Float32}(undef, (r,s,M))

    PsickPIv=Array{Float32}(undef, (r,s,M,V)) # Sum of the infective prodromal individuals, Partially imunized inds stratified by vaccine
    IsickPIv=Array{Float32}(undef, (r,s,M,V)) # Sum of the infective fully infectious individuals, Partially imunized inds stratified by vaccine
    LsickPIv=Array{Float32}(undef, (r,s,M,V)) # Sum of the infective fully infectious individuals, Partially imunized inds stratified by vaccine

    Ql = Array{Float32}(undef, (r))   # symptomatic individuals supposed to be isolated

    for l = 1:r
        Q = 0
        for a = 1:s
            for m = 1:M
              Psum = 0
              Isum = 0
              Lsum = 0
              Isick = 0
              IsickPI = 0
              Lsick = 0
              LsickPI = 0
             # v=V
              for v=1:(V+1) # Individuals waiting to be vaccinated (U)  and vaccinated with pending outcome (V)
                  Psum = Psum + sum(u[INDvecP[l][a][m][v]])
                  Isum = Isum + sum(u[INDvecI[l][a][m][v]])
                  Lsum = Lsum + sum(u[INDvecL[l][a][m][v]])
              end
              # not immunized (NI)
              Psum = Psum + sum(u[INDvecP[l][a][m][2*V+2]])
              Isum = Isum + sum(u[INDvecI[l][a][m][2*V+2]])
              Lsum = Lsum + sum(u[INDvecL[l][a][m][2*V+2]])
              Isick = Isum * fsick[a,m]
              Lsick = Lsum * fsick[a,m]

              for v=1:V # Partially immunized (PI)
                  v1 = v + V+1
                  PErl=sum(u[INDvecP[l][a][m][v1]]) # sum over Erlang states of P
                  PsumPIv[l,a,m,v] = PErl

                  IErl = sum(u[INDvecI[l][a][m][v1]]) # sum over the Erlang states of I
                  IsumPIv[l,a,m,v] = IErl

                  IsickPI = IsickPI + IErl * fsickPI[a][m,v]
                  IsickPIv[l,a,m,v] = IErl * fsickPI[a][m,v]

                  LErl = sum(u[INDvecL[l][a][m][v1]]) # sum over the Erlang states of L
                  LsumPIv[l,a,m,v] = LErl

                  LsickPI = LsickPI + LErl * fsickPI[a][m,v]
                  LsickPIv[l,a,m,v] = LErl * fsickPI[a][m,v]
              end

              PsumU[l,a,m] = Psum
              IsumU[l,a,m] = Isum
              LsumU[l,a,m] = Lsum

              IsickU[l,a,m] = Isick # Symptomatic fully infectious individuals I
              LsickU[l,a,m] = Lsick # Symptomatic late infectious individuals L

              Q = Q + Isick + IsickPI + Lsick + LsickPI # sum of symptomatic individuals going to islolation
            end
        end
        Q = Q * fiso[l]
        is = min(Qmax[l]/Q,1)
        for a = 1:s
            for m = 1 : M

                tmp = (phome + is * (1-phome))* fiso[l] 
                a1 = sum(pPmv[m,:] .* PsumPIv[l,a,m,:])  
                a1 = betaP[a,m] * (PsumU[l,a,m] + a1)
                a2  = sum((1 .- fsickPI[a][m,:] .* tmp) .* pImv[m,:] .* IsumPIv[l,a,m,:]) 
                a2  = betaI[a,m] * (IsumU[l,a,m] - tmp * IsickU[l,a,m] + a2)
                a3  = sum((1 .- fsickPI[a][m,:] .* tmp).* pLmv[m,:] .* LsumPIv[l,a,m,:])
                a3 = betaL[a,m] * (LsumU[l,a,m] - tmp * LsickU[l,a,m] + a3)
                EffInf[(l-1)*s+a,m] = a1 + a2 + a3

            end
        end
    end

    lamex = lamexmat[sum(mutinttimes .<t)]
    lam = zeros(r,s,M)

     for m = 1:M
         lam[:,:,m] = reshape(X * EffInf[:,m],s,r)' * R0new[m]
     end

     if tw1<t<tw2
         t1=t-w_adj[1]
     elseif tw3<t<tw4
        t1=t-w_adj[2]
     else
         t1=t
     end
     lam = lam .* (1+Amp*cos(2*pi*(t1-tR0max)/365)) .+ (lamex)
end

## Model - differential equations
function spatialmodel(du,u,p,t)
    IndSU, IndE1Uam, IndE1U, nEPIL, INDvec, IndLnL, IndRInf = p[1]
    rates,  alpha, vaxrates, vaxtime, fratesS, frates, fE1, fsickdead, delta, betaP, betaI, betaL = p[2]
    N, g, h, R0new, Amp, tR0max = p[3]
    Xmat, XmatEB, XmatNoSch, conttime, tEB, tEBstop,tEB2, tEBstop2 = p[4]
    INDInf = p[5]
    X = p[10]  
    mutint = p[14]
    pPmv,pImv,pLmv = p[15]  
    vrates = vaxrates[sum(vaxtime .< t)+1]
    lamex = p[6]

    lam = lambda(u,INDInf,Qmax,fiso,betaP,betaI,betaL,R0new,Amp,X,tR0max,lamex,t,mutint,pPmv,pImv,pLmv, mutinttimes)/N
    for l = 1:r
        u[IndRInf[s,M,r]+1+l] = 0
        for a = 1:s
            # Susceptibles
            du[IndSU[a,l]] = - u[IndSU[a,l]] * sum(lam[l,a,:]) - u[IndSU[a,l]] * sum(vrates[a,l,:])
            # Incidence unvaccinated (U)
            u[IndRInf[s,M,r]+1+l] += u[IndSU[a,l]]*sum(lam[l,a,:])
            for v = 1:V
                # Sa,l(V,v)
                du[IndSU[a,l]+v] = u[IndSU[a,l]] * vrates[a,l,v] - u[IndSU[a,l]+v] * (alpha[a,v] + sum(lam[l,a,:]))
                # Sa,l(PI,v)
                du[IndSU[a,l]+v+V] = u[IndSU[a,l]+v] * fratesS[v,2] * alpha[a,v]- u[IndSU[a,l]+v+V] * (lam[l,a,:]'g[:,v])
                # Sa,l(Im,v)
                du[IndSU[a,l]+v+2*V] = u[IndSU[a,l]+v] * fratesS[v,1] * alpha[a,v]- u[IndSU[a,l]+v+2*V] * (lam[l,a,:]'h[:,v])
                # Incidence vaccinated (V)
                u[IndRInf[s,M,r]+1+l] += u[IndSU[a,l]+v]*sum(lam[l,a,:]) + u[IndSU[a,l]+v+V]*(lam[l,a,:]'g[:,v]) + u[IndSU[a,l]+v+2*V]*(lam[l,a,:]'h[:,v])
            end

            ## Sa,l(NI)
            du[IndSU[a,l]+3*V+1] = transpose(alpha[a,:] .*fratesS[:,3]) * u[(IndSU[a,l]+1) : (IndSU[a,l]+V)] - sum(lam[l,a,:]) * u[IndSU[a,l]+3*V+1]

            # Incidence NI
            u[IndRInf[s,M,r]+1+l] +=  sum(lam[l,a,:]) * u[IndSU[a,l]+3*V+1]

            for m = 1:M
                ### here E1 compartments
                ## E1uam
                nusum = sum(vrates[a,l,:])

                du[IndE1[l][a][m][1]] = lam[l,a, m] * u[IndSU[a,l]] - (rates[a,m][1] + nusum)  * u[IndE1[l][a][m][1]]
                # E1_a,l(m,V)
                du[IndE1[l][a][m][2:(V+1)]] = lam[l,a, m] .* u[(IndSU[a,l]+1):(IndSU[a,l]+V)]  .+ vrates[a,l,:] .*u[IndE1[l][a][m][1]] - (rates[a,m][1] .+ alpha[a,:]) .* u[IndE1[l][a][m][2:(V+1)]]
                # E1_a,l(m,PI,v)
                du[IndE1[l][a][m][(V+2):(2*V+1)]] = lam[l,a, m]* g[m,:] .* u[(IndSU[a,l]+1+V):(IndSU[a,l]+2*V)] + lam[l,a,m] * h[m,:] .* u[(IndSU[a,l]+1+2*V):(IndSU[a,l]+3*V)] .+ u[IndE1[l][a][m][2:(V+1)]] .* alpha[a,:] .* fE1[m,:] - rates[a,m][1] * u[IndE1[l][a][m][V+2:(2*V+1)]]
                # E1_a,l(m,NI,v)
                du[IndE1[l][a][m][2*V+2]] = lam[l,a, m] .* u[IndSU[a,l]+3*V+1] - u[IndE1[l][a][m][2*V+2]] * rates[a,m][1] + sum(alpha[a,:] .* fE1[m,:] .* u[IndE1[l][a][m][2:(V+1)]])

                ## Infected states
                pick = INDvec[l][a][m]

                ru = rates[a,m] .* u[pick[1][1]]
                du[pick[1][2]] = ru[pick[1][4]] .- ru[pick[1][5]] .- nusum .* u[pick[1][2]] # E2(U) .. L_nL(U)

                ru = rates[a,m] .* u[pick[2*V+2][1]]
                du[pick[2*V+2][2]] = ru[pick[1][4]] - ru[pick[1][5]]  # E2(NI) .. L_nL(NI)  

                RIm = zeros(nEPIL[a,m])
                for v in 2:(V+1)
                    ru = rates[a,m] .* u[pick[v][1]]
                    du[pick[v][2]] = ru[pick[1][4]] - ru[pick[1][5]] + vrates[a,l,v-1] .* u[pick[v][2]]  - alpha[a,v-1] .* u[pick[v][2]] # E2(V) .. L_nL(V)

                    fIm = frates[v-1][a,m][:,1] .* u[pick[v][1]]
                    fPI = frates[v-1][a,m][:,2] .* u[pick[v][1]]
                    fNI = frates[v-1][a,m][:,3] .* u[pick[v][1]]

                    ru = rates[a,m] .* u[pick[V+v][1]]

                    du[pick[V+v][2]] = ru[pick[1][4]] - ru[pick[1][5]] + fPI[pick[1][5]]  # E2(PI) .. L_nL(PI)
                    du[pick[2*V+2][2]] = du[pick[2*V+2][2]] + fNI[pick[1][5]]             # E2(NI) .. L_nL(NI)

                    RIm = RIm .+ fIm
                end
                # R_al (Im,m)
                tmp = 1 .-fsickdeadPI[a][m,:]
                du[IndRInf[a,m,l]] = sum(RIm) + (sum(u[IndLnL[l][a][m][1]]) * (1-fsickdead[a,m])  + u[IndLnL[l][a][m][2]]'tmp) .* delta[a,m]

                # D_al (m)
                du[IndRInf[a,m,l]+1] = ((sum(u[IndLnL[l][a][m][1]]) * fsickdead[a,m])  + u[IndLnL[l][a][m][2]]'fsickdeadPI[a][m,:]) .* delta[a,m]
            end
        end
    end
    du
end

# Callback/event at integer time points for triggered contact reductions
function integer_time_values_condition(u,t,integrator)
    t in integrator.p[7]
end

function contact_reduction_affect!(integrator)
    Sol_t    = integrator.sol.t
    Cur_t    = integrator.t
    Sol_u    = integrator.sol.u
    idx_orig = integrator.p[8]
    xmatr    = integrator.p[4]#[2]
    tschool  = integrator.p[13]
    tEB      = integrator.p[4][5] 
    tEBstop  = integrator.p[4][6]
    tEB2     = integrator.p[4][7]  
    tEBstop2 = integrator.p[4][8]
    if !(tEB <= Cur_t <= tEBstop || tEB2 <= Cur_t <= tEBstop2)
        xmatr1 = xmatr[1]
        integrator.p[10] = xmatr1[sum(xmatr[4] .<= Cur_t)]
    else
        if sum((tschool[:,1] .< Cur_t) .*(tschool[:,2] .> Cur_t)) .> 0
            xmatr1=xmatr[3]
        else
            xmatr1=xmatr[2]
        end
        k1 = sum(Sol_t .<= (Cur_t-7))   
        k2 = sum(Sol_t .<= Cur_t)      
        sel = 1
        for l in 1:integrator.p[9]
            idx_inc = idx_orig+l
            out = zeros((1,k2))
            Incid = 0.

            # Integration of incidence between t-1 and t
            for i in k1:(k2-1)
                out[i] =(Sol_t[i+1] - Sol_t[i])*Sol_u[i+1][idx_inc]
            end
            # Incidence per 100000
            Incid = sum(out)*100000/integrator.p[11]
            # Activating incidence trigger at location l given the incidence
            sel += sum(integrator.p[12] .<= Incid) *5^(l-1)

        end
        integrator.p[10] = xmatr1[sel]
    end
end

###############################################################
###### Solving the IVP
###############################################################

u0 = fill(0.,IndRInf[s,M,r]+1+r)

for l=1:r
    for a = 1:s
        u0[IndSU[a,l]] = max(propvax[l][a]*(PopSize[l,a]-ininf[l][a]),0.)
        u0[IndE1Uam[a,1,l]+nE[a,1]+nP[a,1]] = ininf[l][a]
        u0[IndSU[a,l]+3*V+1] = (1-propvax[l][a])*(PopSize[l,a]-ininf[l][a])
    end
end


# Callbacks for incidence-triggered contact reductions
cb = DiscreteCallback(integer_time_values_condition, contact_reduction_affect!)

# Setting up the ODE problrm for the solver
prob = ODEProblem(spatialmodel,u0,tspan,p,callback=cb)

# Solving the ODE system
sol = solve(prob, VCABM(), tstops=thres)

##########################################################################
###### Post processing solutions to plot the dynamics of all compartments
##########################################################################

# Total number of susceptible individuals
function sumsus(u)
        sumSus=zeros(4)
        for l =1:r  # cummulative infections per location and age group
           for a=1:s
               sumSus=sumSus .+[u[IndSU[a,l]], sum(u[(IndSU[a,l]+1):(IndSU[a,l]+V)]), sum(u[(IndSU[a,l]+V+1):(IndSU[a,l]+3*V)]), u[IndSU[a,l]+3*V+1]]
           end
        end
        sumSus
end

# Total number of infected individuals
function suminf(u)
        sumInf = zeros(M)
        for l = 1:r  # cummulative infections per location and age group
           for a = 1:s
               for m = 1:M
                   for v = 1:2*(V+1)
                       sumInf[m] = sumInf[m] + sum(u[INDvec[l][a][m][v][1]])
                   end
               end
           end
        end
        sumInf
end

# Age-stratified number of infected individuals
function suminfage(u)
    sumInf = zeros(s)
    for l = 1:r  # cummulative infections per location and age group
       for m = 1:M
           for a = 1:s
               for v = 1:2*(V+1)
                   sumInf[a] = sumInf[a] + sum(u[INDvec[l][a][m][v][1]])
               end
           end
       end
    end
    sumInf
end

# Variant-stratified number of infected individuals
function suminfvar(u)
    sumInf = zeros(M)
    for l = 1:r  # cummulative infections per location and age group
       for a = 1:s
           for m = 1:M
               for v = 1:2*(V+1)
                   sumInf[m] = sumInf[m] + sum(u[INDvec[l][a][m][v][1]])
               end
           end
       end
    end
     sumInf
end

# Total number of dead individuals
function sumdead(u)
        sumdead = 0
        for l =1:r  # cummulative infections per location and age group
           for a=1:s
               for m=1:M
                   sumdead += u[IndRInf[a,m,l]+1]
               end
           end
        end
        sumdead
end

# Age-stratified number of dead individuals
function sumdeadage(u)
    sumdead = zeros(s)
    for l = 1:r  # cummulative infections per location and age group
       for a = 1:s
           for m = 1:M
               sumdead[a] = sumdead[a] + u[IndRInf[a,m,l]+1]
           end
       end
    end
    sumdead
end

# Age-stratified number of recovered individuals
function sumrecage(u)
    sumrec=zeros(s)
    for l = 1:r  # cummulative infections per location and age group
       for a = 1:s
           for m = 1:M
               sumrec[a] = sumrec[a] + u[IndRInf[a,m,l]]
           end
       end
    end
    sumrec
end

# Variant-stratified number of dead individuals
function sumdeadvar(u)
    sumdead = zeros(M)
    for l = 1:r  # cummulative infections per location and age group
       for a = 1:s
           for m = 1:M
               sumdead[m] = sumdead[m] + u[IndRInf[a,m,l]+1]
           end
       end
    end
    sumdead
end

# No. of vaccinated individuals
function vacc(u,p,t)
    IndSU, IndE1Uam, IndE1U, nEPIL, INDvec, IndLnL, IndRIm = p[1]
    rates,  alpha, vaxrates, vaxtime, fratesS, frates, fE1, fsickdead, delta, betaP, betaI, betaL = p[2]
    N, g, h, R0new, Amp, tR0max = p[3]
    Xmat, XmatEB, XmatNoSch, conttime, tEB, tEBstop = p[4]
    INDInf = p[5]
    X = p[10] 
    mutint = p[14]
    pPmv,pImv,pLmv = p[15]
           
    vrates = vaxrates[sum(vaxtime .< t)+1]
    du = Array{Union{Missing, Any}}(missing, V*s*r)
    for l = 1:r
        u[IndRIm[s,M,r]+1+l] = 0
        for a = 1:s
            for v = 1:V
                du[V*s*(l-1)+V*(a-1)+v] = u[IndSU[a,l]] * vrates[a,l,v]
            end
            for m = 1:M
                nusum = sum(vrates[a,l,:])
                du[(V*s*(l-1)+V*(a-1)+1):(V*s*(l-1)+V*a)] += vrates[a,l,:] .*u[IndE1[l][a][m][1]]
                pick = INDvec[l][a][m]
                RIm = zeros(nEPIL[a,m])
                for v in 2:(V+1)
                    du[V*s*(l-1)+V*(a-1)+v-1] +=  sum(vrates[a,l,v-1] .* u[pick[v][2]] ) # E2(V) .. L_nL(V)
                end
            end
        end
    end
    du
end

tmx = tmax
ttt = size(sol.u,1)

outInf=fill(0.,(ttt,M))
outInfAgeAll=fill(0.,(ttt,s))
outInfVarAll=fill(0.,(ttt,M))

outDeadAgeAll=fill(0.,(ttt,s))
outRecAgeAll=fill(0.,(ttt,s))
outDeadVarAll=fill(0.,(ttt,M))

outInfAll=fill(0.,(ttt,1))
outSus=fill(0.,(ttt,4))
outDead=fill(0.,ttt)
outInc=fill(0.,ttt)
outIncl1=fill(0.,ttt)
outIncl2=fill(0.,ttt)

outVacc=fill(0.,ttt,V*r*s)
outVaccAll=fill(0.,ttt,1)
outVaccAge=fill(0.,ttt,s)
outVac=fill(0.,tmx,1)
outVacAge=fill(0.,tmx,s)
outVaccCum=fill(0.,ttt,V*r*s)

for k=1:size(sol.u,1)
    outSus[k,:]=sumsus(sol.u[k])
    outInf[k,:]=suminf(sol.u[k])
    outInfAgeAll[k,:]=suminfage(sol.u[k])
    outDeadAgeAll[k,:]=sumdeadage(sol.u[k])
    outRecAgeAll[k,:]=sumrecage(sol.u[k])
    outDeadVarAll[k,:]=sumdeadvar(sol.u[k])
    outInfVarAll[k,:]=suminfvar(sol.u[k])
    outInfAll[k]=sum(outInf[k,:])
    outInc[k]=sum(sol.u[k][(IndRInf[s,M,r]+1+1):(IndRInf[s,M,r]+1+r)])
    outIncl1[k]=sol.u[k][(IndRInf[s,M,r]+1+1)]
    outDead[k]=sumdead(sol.u[k])
end

tmx=tmax
outInc1=fill(0.,tmx,1)
outInc2=fill(0.,tmx,1)

outInc1l1=fill(0.,tmx,1)
outInc2l1=fill(0.,tmx,1)

outInc1l2=fill(0.,tmx,1)
outInc2l2=fill(0.,tmx,1)

outDead1=fill(0.,tmx,1)
outt=fill(0.,tmx,1)

outInf=fill(0.,tmx,1)
outInfAge=fill(0.,tmx,s)
outInfVar=fill(0.,tmx,M)
outdeadAge=fill(0.,tmx,s)
outRecAge=fill(0.,tmx,s)
outdeadVar=fill(0.,tmx,M)

# Numerical integration of incidence
k1 = 1
for t in 1:(tmx)
        k2 = sum(sol.t .<=t)
        inc=0.
        inc += sum((sol.t[(k1+1):k2] .-sol.t[k1:(k2-1)]) .* outInc[(k1+1):k2])

        inc1=0.
        inc1 += sum((sol.t[(k1+1):k2] .-sol.t[k1:(k2-1)]) .* outIncl1[(k1+1):k2])

        inc2=0.
        inc2 += sum((sol.t[(k1+1):k2] .-sol.t[k1:(k2-1)]) .* outIncl2[(k1+1):k2])

        outInc1l1[t,1] = inc1

        outInc1l2[t,1] = inc2

        outInc1[t,1] = inc
        outt[t,1] = sol.t[k2]
        k1=k2
        outDead1[t] = outDead[k2]
end

# 7-days average of new infetions incidence
for t in 1:tmx
        if t>7
            outInc2[t,1] = sum(outInc1[(t-7):t,1])/7
            outInc2l1[t,1] = sum(outInc1l1[(t-7):t,1])/7
            outInc2l2[t,1] = sum(outInc1l2[(t-7):t,1])/7
            else
            outInc2[t,1] = sum(outInc1[1:t,1])/t
            outInc2l1[t,1] = sum(outInc1l1[1:t,1])/t
            outInc2l2[t,1] = sum(outInc1l2[1:t,1])/t
        end
end

for t in 1:tmx
        k = sum(sol.t .<= t)
        outInf[t] = outInfAll[k]
end

for t in 1:tmx
    k = sum(sol.t .<= t)
    outInfAge[t,:] = outInfAgeAll[k,:]
    outdeadAge[t,:] = outDeadAgeAll[k,:]
    outRecAge[t,:] = outRecAgeAll[k,:]
    outInfVar[t,:] = outInfVarAll[k,:]
    outdeadVar[t,:] = outDeadVarAll[k,:]
end

ttt = size(sol.u,1)
deltat = sol.t[2:ttt]-sol.t[1:(ttt-1)]

for k = 1:size(sol.u,1)
    outVacc[k,:] = vacc(sol.u[k],p,sol.t[k])
end

for a = 1:s
    for v=1:V
        outVaccCum[:,(V-1)*a+v] = cumsum( vcat(0,deltat .* outVacc[2:(ttt),((V-1)*a+v)]))
    end
    outVaccAge[:,a] =  sum(outVaccCum[:,((a-1)*V+1):a*V], dims=2)
end

outVaccAll = sum(outVaccAge,dims=2)
for t in 1:tmx
    k = sum(sol.t .<= t)
    outVac[t] = outVaccAll[k]
end

for t in 1:tmx
    k = sum(sol.t .<= t)
    outVacAge[t,:] = outVaccAge[k,:]
end

############################################################################
###### Saving the outputs as a dataframe 
######
###### outt:       Time points
###### outInc2l1:  7-day average incidence in location l1
###### outdead:    Number of dead individuals
###### outVac:     Number of vaccinated individuals
###### outInfAge:  Number of infections per age group
###### outInfVar:  Number of infections per viral virant
###### outdeadAge: Number of dead individuals per age
###### outdeadVar: Number of dead individuals per viral variant
###### outVacAge:  Number of vaccinated individuals per age group
############################################################################

out1 = hcat(outt, outInc2, outInf, outDead1, outVac, outInfAge, outInfVar, outdeadAge, outdeadVar, outVacAge)
out = DataFrame(Tables.table(out1))

# Enter path where the dataframe should be saved (replace "foo/" by the path to the location you want to save)
path = "foo/"

# Enter name of the saved file (replace "test" by the name you want your file to have)
namefile = "test"

# Save the file in the .txt format
CSV.write(path*namefile*".txt", out, header=false, append=false)
