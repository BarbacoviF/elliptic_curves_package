# Curve parameters ------------------------------------------------------------------------------------------------------
# Taken from [https://github.com/arkworks-rs/curves/tree/master/mnt4_753/src/curves]

# Seed
u = -0x15474b1d641a3fd86dcbcee5dcda7fe51852c8cbe26e600733b714aa43c31a66b0344c4e2c428b07a7713041ba18000

# Signed base two decomposition of abs(u) - MSB to LSB (minus_exp_u is LSB to MSB)
minus_exp_miller_loop = [
        1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0, -1, 0, 1, 0, 1, 0, -1, 0, -1, 0, 0, 1, 0, 0, 0, -1, 0,
        -1, 0, -1, 0, 0, 1, 0, 0, 0, 0, 1, 0, -1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, -1,
        0, 0, 0, 1, 0, 0, -1, 0, 0, -1, 0, 0, -1, 0, 1, 0, -1, 0, 0, 0, -1, 0, 1, 0, 0, 0, -1, 0,
        0, -1, 0, 1, 0, -1, 0, 0, 0, -1, 0, 0, -1, 0, 1, 0, 0, -1, 0, -1, 0, 1, 0, 1, 0, 0, 0, 0,
        0, 0, 0, 0, 0, -1, 0, 0, 1, 0, 1, 0, 0, 1, 0, -1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, -1, 0, -1,
        0, 0, 1, 0, 0, 1, 0, -1, 0, 1, 0, -1, 0, 0, 0, 0, -1, 0, 0, 0, 1, 0, 1, 0, 0, -1, 0, 0, -1,
        0, 1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, -1, 0, 1, 0, -1, 0, 1, 0, 0, 0, -1, 0, 0,
        -1, 0, 0, -1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0,
        -1, 0, 0, 0, 1, 0, -1, 0, 0, 1, 0, -1, 0, 1, 0, 1, 0, -1, 0, 1, 0, 0, -1, 0, -1, 0, -1, 0,
        0, 0, 0, 0, 1, 0, -1, 0, 1, 0, 0, 0, 1, 0, 1, 0, -1, 0, 0, 0, 1, 0, 1, 0, 0, -1, 0, 0, 1,
        0, -1, 0, -1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, -1, 0, -1, 0, 0, 0, 0, 1, 0, 0,
        0, -1, 0, 1, 0, 1, 0, 0, 0, -1, 0, 0, -1, 0, 0, 0, 1, 0, 1, 0, -1, 0, 0, 0, 0, 0, 1, 0, 0,
        0, 0, 1, 0, 0, -1, 0, 0, -1, 0, 1, 0, 0, 0, 1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0,
    ][::-1]
exp_miller_loop = [-el for el in minus_exp_miller_loop]

# Modulus
q = u**2 + u + 1

# r-torsion = q - t + 1
r = u**2 + 1

# Curve coefficients
a = 2
b = 28798803903456388891410036793299405764940372360099938340752576406393880372126970068421383312482853541572780087363938442377933706865252053507077543420534380486492786626556269083255657125025963825610840222568694137138741554679540

# Non-residues
NON_RESIDUE_FQ = [13]   # List serialisation
NON_RESIDUE_FQ2 = [0,1] # List serialisation

# Cofactors
h1 = 1
h2 = 41898490967918953402344214791240637128170709919953949071783502921025352812571106773058893763790338921418070971888049094905534395567574915333486969589229856772141392370549616644545554517640527237829320384324374366385444967219201

# Generators
g1_X = 7790163481385331313124631546957228376128961350185262705123068027727518350362064426002432450801002268747950550964579198552865939244360469674540925037890082678099826733417900510086646711680891516503232107232083181010099241949569
g1_Y = 6913648190367314284606685101150155872986263667483624713540251048208073654617802840433842931301128643140890502238233930290161632176167186761333725658542781350626799660920481723757654531036893265359076440986158843531053720994648
g2_X0 = 29483965110843144675703364744708836524643960105538608078862508397502447349913068434941060515343254862580437318493682762113105361632548148204806052114008731372757389645383891982211245013965175213456066452587869519098351487925167
g2_X1 = 19706011319630172391076079624799753948158506771222147486237995321925443331396169656568431378974558350664383559981183980668976846806019030432389169137953988990802000581078994008283967768348275973921598166274857631001635633631000
g2_Y0 = 39940152670760519653940320314827327941993141403708338666925204282084477074754642625849927569427860786384998614863651207257467076192649385174108085803168743803491780568503369317093191101779534035377266300185099318717465441820654
g2_Y1 = 17608637424964395737041291373756657139607306440193731804102457011726690702169238966996114255971643893157857311132388792357391583164125870757541009035041469463366528798593952884745987697403056488744603829437448927398468360797245

# -------------------------------------------------------------------------------------------------------------------

