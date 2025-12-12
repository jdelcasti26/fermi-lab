% H2 model validated PYSCF
% 4 qubits
HH = pauliString.empty;
HH(1, 1) = pauliString(o.PDH, [-0.812170607248712], {'IIII'});
HH(1, 2) = pauliString(o.PDH, [-0.0453026155037992], {'YYXX'});
HH(1, 3) = pauliString(o.PDH, [0.0453026155037992], {'XYYX'});
HH(1, 4) = pauliString(o.PDH, [0.0453026155037992], {'YXXY'});
HH(1, 5) = pauliString(o.PDH, [-0.0453026155037992], {'XXYY'});
HH(1, 6) = pauliString(o.PDH, [0.171412826447769], {'IIIZ'});
HH(1, 7) = pauliString(o.PDH, [0.168688981703612], {'IIZZ'});
HH(1, 8) = pauliString(o.PDH, [0.120625234833904], {'IZIZ'});
HH(1, 9) = pauliString(o.PDH, [0.165927850337703], {'ZIIZ'});
HH(1, 10) = pauliString(o.PDH, [0.171412826447769], {'IIZI'});
HH(1, 11) = pauliString(o.PDH, [0.165927850337703], {'IZZI'});
HH(1, 12) = pauliString(o.PDH, [0.120625234833904], {'ZIZI'});
HH(1, 13) = pauliString(o.PDH, [-0.223431536908134], {'IZII'});
HH(1, 14) = pauliString(o.PDH, [0.174412876122615], {'ZZII'});
HH(1, 15) = pauliString(o.PDH, [-0.223431536908134], {'ZIII'});


