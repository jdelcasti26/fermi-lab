% Ising model Open boundaries on 6 qubits
HH = pauliString.empty;
HH(1, 1) = pauliString(o.PDH, [-1], {'IIIIIX'});
HH(1, 2) = pauliString(o.PDH, [-1], {'IIIIXI'});
HH(1, 3) = pauliString(o.PDH, [-1], {'IIIXII'});
HH(1, 4) = pauliString(o.PDH, [-1], {'IIXIII'});
HH(1, 5) = pauliString(o.PDH, [-1], {'IXIIII'});
HH(1, 6) = pauliString(o.PDH, [-1], {'XIIIII'});
HH(1, 7) = pauliString(o.PDH, [-1], {'IIIIZZ'});
HH(1, 8) = pauliString(o.PDH, [-1], {'IIIZZI'});
HH(1, 9) = pauliString(o.PDH, [-1/2], {'IIZZII'});
HH(1, 10) = pauliString(o.PDH, [-1/2], {'IIZZII'});
HH(1, 11) = pauliString(o.PDH, [-1], {'IZZIII'});
HH(1, 12) = pauliString(o.PDH, [-1], {'ZZIIII'});
