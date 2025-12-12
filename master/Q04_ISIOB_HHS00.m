% Ising model Open boundaries on 4 qubits
HH = pauliString.empty;
HH(1, 1) = pauliString(o.PDH, [-1], {'IIIX'});
HH(1, 2) = pauliString(o.PDH, [-1], {'IIXI'});
HH(1, 3) = pauliString(o.PDH, [-1], {'IXII'});
HH(1, 4) = pauliString(o.PDH, [-1], {'XIII'});
HH(1, 5) = pauliString(o.PDH, [-1], {'IIZZ'});
HH(1, 6) = pauliString(o.PDH, [-1], {'IZZI'});
HH(1, 7) = pauliString(o.PDH, [-1], {'ZZII'});
