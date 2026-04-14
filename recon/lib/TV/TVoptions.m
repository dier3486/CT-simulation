function TVopts = TVoptions(TVopts)

arguments
    TVopts.DIM (1,1) {mustBeMember(TVopts.DIM, ["1D", "2D", "3D"])} = "3D"
    TVopts.MaxIter (1,1) {mustBeNumeric} = 100
    TVopts.tol (1,1) {mustBeNumeric} = 0.01
    TVopts.Crange (1,2) {mustBeNumeric} = [-inf inf]
    TVopts.Balance (1,:) {mustBeNumeric} = [1 1 1]
    TVopts.Pursmit (1,1) {mustBeMember(TVopts.Pursmit, ["I", "tanh", "log"])} = "I"
    TVopts.PursmitCoeff (1,:) {mustBeNumeric} = []
    TVopts.Kinetic (1,1) {mustBeMember(TVopts.Kinetic, ["0", "Laplace"])} = "0"
    TVopts.KineticCoeff (1,:) {mustBeNumeric} = []
    TVopts.LeftBoundary (1,1) {mustBeNumericOrLogical} = false
    TVopts.LeftDim (1,:) {mustBeNumeric} = 1
    TVopts.BoudaryCondX (1,1) ...
        {mustBeMember(TVopts.BoudaryCondX, ["fix", "clamp", "linextrap", "mirror", "oddextrap"])} = "clamp"
    TVopts.BoudaryCondY (1,1) ...
        {mustBeMember(TVopts.BoudaryCondY, ["fix","clamp", "linextrap", "mirror", "oddextrap"])} = "clamp"
    TVopts.BoudaryCondZ (1,1) ...
        {mustBeMember(TVopts.BoudaryCondZ, ["fix","clamp", "linextrap", "mirror", "oddextrap"])} = "clamp"
end

% normalize the Balance
Ndim = str2double(regexp(TVopts.DIM, '^\d', 'match'));
balance = TVopts.Balance;
if length(balance) < Ndim
    balance(end+1 : Ndim) = 1;
else
    balance = balance(1 : Ndim);
end
balance = balance.*(Ndim/sum(balance));
if length(balance) < 3
    balance(end+1 : 3) = 0;
end
TVopts.Balance = balance;

end