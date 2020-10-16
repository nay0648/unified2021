function Demix = heig2(diagload, C1, C0)
%
% solve the generalized eigenvalue decomposition problem for 2x2
% Hermitian matrices for all frequency bins: C1*w = lambda*C0*w,
% the theoretical better solution
% diagload:             the diagonal loading to prevent singular
% C1:                   the matrix C1, 2 x 2
% C0:                   the matrix C0, 2 x 2
% Demix                 the demixing matrix, rows are the eigenvectors
%                       the scaling problem is solved
%

d = real(C1(1, 1)) + diagload;
er = real(C1(1, 2));
ei = imag(C1(1, 2));
f = real(C1(2, 2)) + diagload;

a = real(C0(1, 1)) + diagload;
br = real(C0(1, 2));
bi = imag(C0(1, 2));
c = real(C0(2, 2)) + diagload;

%
% calculate eigenvalues
%
x = a * c - (br * br + bi * bi);
y = 2.0 * (br * er + bi * ei) - d * c - a * f;
z = d * f - (er * er + ei * ei);

delta = sqrt(y * y - 4.0 * x * z);
lambda0 = (-y + delta) / (2.0 * x);
lambda1 = (-y - delta) / (2.0 * x);

%
% calculate the eigenvectors, as the row vectors of B
%
b00 = f - lambda0 * c;
b01 = (-er + lambda0 * br) + (-ei + lambda0 * bi) * 1i;
b10 = (-er + lambda1 * br) + (ei - lambda1 * bi) * 1i;
b11 = d - lambda1 * a;

%
% solve the scaling ambiguity
%

%
% calculate A = B^-1
%
detbr = (real(b00) * real(b11) - imag(b00) * imag(b11)) - (real(b01) * real(b10) - imag(b01) * imag(b10));
detbi = (real(b00) * imag(b11) + imag(b00) * real(b11)) - (real(b01) * imag(b10) + imag(b01) * real(b10));

detbsq = detbr * detbr + detbi * detbi;
detb = (detbr / detbsq) - (detbi / detbsq) * 1i;

a00 = detb * b11;
a11 = detb * b00;
a01 = -detb * b01;
a10 = -detb * b10;

if abs(a00) >= abs(a10)
    Demix(1, 1) = a00 * b00;
    Demix(1, 2) = a00 * b01;
else
    Demix(1, 1) = a10 * b00;
    Demix(1, 2) = a10 * b01;
end

if abs(a11) >= abs(a01)
    Demix(2, 1) = a11 * b10;
    Demix(2, 2) = a11 * b11;
else
    Demix(2, 1) = a01 * b10;
    Demix(2, 2) = a01 * b11;
end

end
