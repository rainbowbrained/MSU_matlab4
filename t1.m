% Œ 
% Œ 
clear, clc;
mex -setup C++
mex cubesolve.cpp

Acompl = [11+1.0*1i 2+1.0*1i];
Bcompl = [10+1.0*1i 1+1.0*1i];
Ccompl = [1+1.0*1i 1+1.0*1i];

[x1] = cubesolve(Acompl, Bcompl, Ccompl)
% A(x-x1)(x-x2) = A(x^2 - (x1+x2)x +x1x2)
Acompl.*[x1; x2].^3 + Bcompl.*[x1; x2] + Ccompl
A = [0 0 0];
B = [0 3 2];
C = [1 2 0];

[x1, x2, D] = cubesolve(A, B, C)
% A(x-x1)(x-x2) = A(x^2 - (x1+x2)x +x1x2)
A.*[x1; x2].^2 + B.*[x1; x2] + C
pause

clear, clc;
mex -setup C++
mex mexfunc.cpp

Acompl = [11+1.0*1i 2+1.0*1i];
Bcompl = [10+1.0*1i 1+1.0*1i];
Ccompl = [1+1.0*1i 1+1.0*1i];

[x1, x2, D] = mexfunc(Acompl, Bcompl, Ccompl)
% A(x-x1)(x-x2) = A(x^2 - (x1+x2)x +x1x2)
Acompl.*[x1; x2].^2 + Bcompl.*[x1; x2] + Ccompl
A = [0 0 0];
B = [0 3 2];
C = [1 2 0];

[x1, x2, D] = mexfunc(A, B, C)
% A(x-x1)(x-x2) = A(x^2 - (x1+x2)x +x1x2)
A.*[x1; x2].^2 + B.*[x1; x2] + C