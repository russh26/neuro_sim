function CreateInxCoefficients
Params = load_csv('sFile','MileusnicParams.csv');

%% Bag1

ksr(1) = Params.Bag1(1);
kpr(1) = Params.Bag1(2);
m(1) = Params.Bag1(3);
Beta0(1) = Params.Bag1(4);
Beta1(1) = Params.Bag1(5);
Beta2(1) = Params.Bag1(6);
g1(1) = Params.Bag1(7);
g2(1) = Params.Bag1(8);
cl(1) = Params.Bag1(9);
cs(1) = Params.Bag1(10);
X(1) = Params.Bag1(11);
lsr(1) = Params.Bag1(12);
lpr(1) = Params.Bag1(13);
G(1) = Params.Bag1(14);
a(1) = Params.Bag1(15);
R(1) = Params.Bag1(16);
l0sr(1) = Params.Bag1(17);
l0pr(1) = Params.Bag1(18);
Lsec(1) = Params.Bag1(19);
tau(1) = Params.Bag1(20);
freq(1) = Params.Bag1(21);
p(1) = Params.Bag1(22);

%% Bag2

ksr(2) = Params.Bag2(1);
kpr(2) = Params.Bag2(2);
m(2) = Params.Bag2(3);
Beta0(2) = Params.Bag2(4);
Beta1(2) = Params.Bag2(5);
Beta2(2) = Params.Bag2(6);
g1(2) = Params.Bag2(7);
g2(2) = Params.Bag2(8);
cl(2) = Params.Bag2(9);
cs(2) = Params.Bag2(10);
X(2) = Params.Bag2(11);
lsr(2) = Params.Bag2(12);
lpr(2) = Params.Bag2(13);
G(2) = Params.Bag2(14);
a(2) = Params.Bag2(15);
R(2) = Params.Bag2(16);
l0sr(2) = Params.Bag2(17);
l0pr(2) = Params.Bag2(18);
Lsec(2) = Params.Bag2(19);
tau(2) = Params.Bag2(20);
freq(2) = Params.Bag2(21);
p(2) = Params.Bag2(22);

%% Chain

ksr(3) = Params.Chain(1);
kpr(3) = Params.Chain(2);
m(3) = Params.Chain(3);
Beta0(3) = Params.Chain(4);
Beta1(3) = Params.Chain(5);
Beta2(3) = Params.Chain(6);
g1(3) = Params.Chain(7);
g2(3) = Params.Chain(8);
cl(3) = Params.Chain(9);
cs(3) = Params.Chain(10);
X(3) = Params.Chain(11);
lsr(3) = Params.Chain(12);
lpr(3) = Params.Chain(13);
G(3) = Params.Chain(14);
a(3) = Params.Chain(15);
R(3) = Params.Chain(16);
l0sr(3) = Params.Chain(17);
l0pr(3) = Params.Chain(18);
Lsec(3) = Params.Chain(19);
tau(3) = Params.Chain(20);
freq(3) = Params.Chain(21);
p(3) = Params.Chain(22);

clear Params
save('Coefficients');