% PROTOTYPE TEST
% ARRAY OF PAULI STRINGS WITH UCCSD ANTIHERMITIAN OPERATORS
AZ = pauliString.empty;
AZ(1, 1) = pauliString(o.PDH, [-0.5j; +0.5j], {'IYZX';'IXZY'});
AZ(1, 2) = pauliString(o.PDH, [-0.5j; +0.5j], {'YZXI';'XZYI'});
AZ(1, 3) = pauliString(o.PDH, [-0.125j;-0.125j;0.125j;-0.125j;0.125j;-0.125j;0.125j;0.125j],...
                              {'YXXX';'XYXX';'XXYX';'YYYX';'XXXY';'YYXY';'YXYY';'XYYY'});

%Individual T - T^ singles
%-1 [0^ 2] +
%1 [2^ 0]
%-0.5j [X0 Z1 Y2] +
%0.5j [Y0 Z1 X2]
%Individual T - T^ singles
%-1 [1^ 3] +
%1 [3^ 1]
%-0.5j [X1 Z2 Y3] +
%0.5j [Y1 Z2 X3]
%Individual T - T^ doubles
%1 [1^ 0^ 3 2] +
%-1 [3^ 2^ 1 0]
%-0.125j [X0 X1 X2 Y3] +
%-0.125j [X0 X1 Y2 X3] +
%0.125j [X0 Y1 X2 X3] +
%-0.125j [X0 Y1 Y2 Y3] +
%0.125j [Y0 X1 X2 X3] +
%-0.125j [Y0 X1 Y2 Y3] +
%0.125j [Y0 Y1 X2 Y3] +
%0.125j [Y0 Y1 Y2 X3]