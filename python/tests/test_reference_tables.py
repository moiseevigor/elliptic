"""Cross-language reference tests: exact numerical values taken from the
MATLAB/Octave test suite (matlab/tests/test*.m).  Every constant here was
produced by the Octave implementation and verified to < 1e-12 against
independent sources (DLMF tables, scipy.special).

Any disagreement between Python and MATLAB at these tolerances is a regression.
"""
import numpy as np
import pytest
import elliptic


# ---------------------------------------------------------------------------
# elliptic12 — 5×7 meshgrid (phi in 0:15:90 deg, alpha in 0:20:90 deg)
# Taken verbatim from matlab/tests/testElliptic12.m
# ---------------------------------------------------------------------------

_F12_PHI = np.deg2rad(np.array([0, 15, 30, 45, 60, 75, 90], dtype=float))
_F12_M   = np.sin(np.deg2rad(np.array([[0], [20], [40], [60], [80]], dtype=float))) ** 2

_EXPECTED_F = np.array([
    [0, 0.261799387799149, 0.523598775598299, 0.785398163397448, 1.047197551196598, 1.308996938995747, 1.570796326794897],
    [0, 0.262145681692449, 0.526283990562203, 0.793981429961732, 1.065968913708522, 1.341839009622797, 1.620025899124204],
    [0, 0.263033690353369, 0.533427451037688, 0.818147652199543, 1.122556669683242, 1.447669376453414, 1.786769134885021],
    [0, 0.264063548276829, 0.542229109803553, 0.851223749071185, 1.212596615254979, 1.649178665655556, 2.156515647499643],
    [0, 0.264747663195422, 0.548425344542772, 0.877408330405700, 1.301353213761152, 1.946822305295344, 3.153385251887838],
])

_EXPECTED_E = np.array([
    [0, 0.261799387799149, 0.523598775598299, 0.785398163397448, 1.047197551196598, 1.308996938995747, 1.570796326794897],
    [0, 0.261453912906691, 0.520937696463332, 0.776974018512720, 1.028972213953050, 1.277421529463950, 1.523799205259775],
    [0, 0.260575450957998, 0.514088617513829, 0.754888085407355, 0.980134299660981, 1.191010358808463, 1.393140248523812],
    [0, 0.259569955098381, 0.506092072465726, 0.728224155457347, 0.918393294316326, 1.075856687976766, 1.211056027568461],
    [0, 0.258909826877359, 0.500742319367686, 0.709723805114461, 0.872755203912916, 0.981407813910667, 1.040114395706011],
])

_EXPECTED_Z = np.array([
    [0, 0,                  0,                  0,                  0,                  0,                  0],
    [0, 0.014879224414869,  0.025914050857933,  0.030153691368092,  0.026319982023431,  0.015285400927267,  0.000000000000000],
    [0, 0.055488619316143,  0.098176322411480,  0.116980030437629,  0.104879154913245,  0.062265499992439,  0.000000000000000],
    [0, 0.111277151300479,  0.201587056463619,  0.250193934198656,  0.237423303852160,  0.149710955672358,  0.000000000000000],
    [0, 0.171585306172434,  0.319849389939710,  0.420318939400386,  0.443516115310938,  0.339266830849054,  0.000000000000000],
])


class TestElliptic12Grid:
    def test_F_grid(self):
        F, _, _ = elliptic.elliptic12(_F12_PHI, _F12_M)
        np.testing.assert_allclose(F, _EXPECTED_F, atol=1e-12, err_msg="F grid mismatch vs MATLAB")

    def test_E_grid(self):
        _, E, _ = elliptic.elliptic12(_F12_PHI, _F12_M)
        np.testing.assert_allclose(E, _EXPECTED_E, atol=1e-12, err_msg="E grid mismatch vs MATLAB")

    def test_Z_grid(self):
        _, _, Z = elliptic.elliptic12(_F12_PHI, _F12_M)
        np.testing.assert_allclose(Z, _EXPECTED_Z, atol=1e-12, err_msg="Z grid mismatch vs MATLAB")

    def test_large_argument(self):
        """testElliptic12.m: elliptic12(1000*pi/e, 0.5)"""
        F, E, Z = elliptic.elliptic12(1000 * np.pi / np.e, 0.5)
        assert abs(float(F) - 1364.215673994739) < 1e-8
        assert abs(float(E) - 993.6995958059659) < 1e-8
        assert abs(float(Z) - (-9.508521098575250e-02)) < 1e-8


# ---------------------------------------------------------------------------
# elliptic3 — spot value and 3-D meshgrid
# Taken from matlab/tests/testElliptic3.m
# ---------------------------------------------------------------------------

_EXPECTED_PI = np.array([
    [0, 0, 0],
    [0, 0, 0],
    [0, 0, 0],
    [0, 0, 0],
    [0.436332312998582, 0.447481674284817, 0.459719055660277],
    [0.438747923080541, 0.450008007383737, 0.462369007025916],
    [0.444551505567130, 0.456080142788474, 0.468741143521980],
    [0.449816440114804, 0.461591818685859, 0.474528668429618],
    [0.872664625997165, 0.962368261625533, 1.094942500776336],
    [0.890543879394498, 0.983487876637982, 1.121161300350779],
    [0.940075683068702, 1.042322814335382, 1.194760969067305],
    [0.997105354291443, 1.110665166734768, 1.281307083403909],
    [1.308996938995747, 1.597941344935833, 2.305386400623875],
    [1.360834669550633, 1.668946178349462, 2.428736708263183],
    [1.534546187668256, 1.911003257311173, 2.863258220259662],
    [1.871453962422796, 2.397752070491359, 3.803705851444642],
])


class TestElliptic3Grid:
    def test_spot_value(self):
        """pi/4, m=0.5, n=0.5 — exact from MATLAB testElliptic3.m"""
        Pi = elliptic.elliptic3(np.pi / 4, 0.5, 0.5)
        assert abs(float(Pi) - 0.919022739165697) < 1e-10

    def test_3d_grid(self):
        phi_deg = np.array([0, 25, 50, 75], dtype=float)
        alpha_deg = np.array([0, 25, 50, 75], dtype=float)
        c_vals = np.array([0.0, 0.4, 0.8], dtype=float)

        phi_r = np.deg2rad(phi_deg)
        m_vals = np.sin(np.deg2rad(alpha_deg)) ** 2

        rows = []
        for c in c_vals:
            for phi in phi_r:
                for m in m_vals:
                    rows.append(float(elliptic.elliptic3(phi, m, c)))

        # reshape to match MATLAB column-major order (phi varies fastest, then alpha, then c)
        Pi = np.array(rows).reshape(3, 4, 4).transpose(1, 2, 0).reshape(-1, 3)
        np.testing.assert_allclose(Pi, _EXPECTED_PI, atol=1e-12, err_msg="elliptic3 3D grid mismatch vs MATLAB")


# ---------------------------------------------------------------------------
# ellipticBDJ — B / D / J reference tables (3 m × 4 phi)
# Taken verbatim from matlab/tests/testEllipticBDJ.m (Group 4)
# ---------------------------------------------------------------------------

_BDJ_PHI = np.array([np.pi/6, np.pi/4, np.pi/3, np.pi/2])
_BDJ_M   = np.array([[0.2], [0.5], [0.8]])
_BDJ_PHI_GRID, _BDJ_M_GRID = np.meshgrid(_BDJ_PHI, _BDJ_M)

_EXPECTED_B = np.array([
    [0.4822319097177440, 0.6529667634636220, 0.7569756041614870, 0.8066808960371530],
    [0.4884759118954500, 0.6703551321060780, 0.7874738572402080, 0.8472130847939790],
    [0.4952048001091030, 0.6909260550736680, 0.8281306941693560, 0.9088110737045840],
])

_EXPECTED_D = np.array([
    [4.600298959653437e-02, 1.474128859214848e-01, 3.234022020908615e-01, 8.529427025733753e-01],
    [4.714682090995279e-02, 1.556627441431676e-01, 3.549552008055692e-01, 1.006861592507392e+00],
    [4.839998016351204e-02, 1.658742182004513e-01, 4.016987480555824e-01, 1.348394253116269e+00],
])

_EXPECTED_J = np.array([
    [4.823487054607033e-02, 1.630222842878820e-01, 3.807339316051869e-01, 1.112925198627194],
    [4.944500507300179e-02, 1.723183252575266e-01, 4.189753046197145e-01, 1.321007148808583],
    [5.077140114947837e-02, 1.838496301039670e-01, 4.759688499155832e-01, 1.788630886401574],
])


class TestEllipticBDJTable:
    def test_B_table(self):
        B, _, _ = elliptic.ellipticBDJ(_BDJ_PHI_GRID, _BDJ_M_GRID, 0.3)
        np.testing.assert_allclose(np.asarray(B), _EXPECTED_B, atol=1e-12, err_msg="B table mismatch vs MATLAB")

    def test_D_table(self):
        _, D, _ = elliptic.ellipticBDJ(_BDJ_PHI_GRID, _BDJ_M_GRID, 0.3)
        np.testing.assert_allclose(np.asarray(D), _EXPECTED_D, atol=1e-12, err_msg="D table mismatch vs MATLAB")

    def test_J_table(self):
        _, _, J = elliptic.ellipticBDJ(_BDJ_PHI_GRID, _BDJ_M_GRID, 0.3)
        np.testing.assert_allclose(np.asarray(J), _EXPECTED_J, atol=1e-12, err_msg="J table mismatch vs MATLAB")


# ---------------------------------------------------------------------------
# jacobiEDJ — Eu / Du / Ju tables at u = K/4, K/2, 3K/4, K
# Taken from matlab/tests/testEllipticBDJ.m (Group 4, jacobiEDJ section)
# ---------------------------------------------------------------------------

_JEDJ_M = [0.2, 0.5, 0.8]

_EXPECTED_EU = np.array([
    [0.4103348919732966, 0.7973039335479685, 1.153380434524441, 1.489035058095853],
    [0.4479288836359817, 0.8217685499305638, 1.110593848820260, 1.350643881047676],
    [0.5213087503809286, 0.8656381644139404, 1.055687219771291, 1.178489924327838],
])

_EXPECTED_DU = np.array([
    [0.02285503839667675, 0.1625393287864774, 0.4566863221672750, 0.8529427025733753],
    [0.03117957137872241, 0.2105375774402440, 0.5599243183115373, 1.006861592507392],
    [0.05374072665535606, 0.3287056237456081, 0.7965209691804365, 1.348394253116269],
])

_EXPECTED_JU = np.array([
    [0.02354858117828542, 0.1809070499958038, 0.5570048420102275, 1.112925198627194],
    [0.03233444406148652, 0.2379084919244906, 0.6917734853508237, 1.321007148808583],
    [0.05656046648375568, 0.3828083834319858, 1.007421432719784,  1.788630886401574],
])


def _K(m):
    from scipy.special import ellipk
    return float(ellipk(m))


class TestJacobiEDJTable:
    @pytest.mark.parametrize("mi,m", enumerate(_JEDJ_M))
    def test_Eu_row(self, mi, m):
        K = _K(m)
        u = np.array([0.25*K, 0.5*K, 0.75*K, K])
        Eu, _, _ = elliptic.jacobiEDJ(u, m)
        np.testing.assert_allclose(np.asarray(Eu), _EXPECTED_EU[mi], atol=1e-12,
                                   err_msg=f"Eu table row m={m} vs MATLAB")

    @pytest.mark.parametrize("mi,m", enumerate(_JEDJ_M))
    def test_Du_row(self, mi, m):
        K = _K(m)
        u = np.array([0.25*K, 0.5*K, 0.75*K, K])
        _, Du, _ = elliptic.jacobiEDJ(u, m)
        np.testing.assert_allclose(np.asarray(Du), _EXPECTED_DU[mi], atol=1e-12,
                                   err_msg=f"Du table row m={m} vs MATLAB")

    @pytest.mark.parametrize("mi,m", enumerate(_JEDJ_M))
    def test_Ju_row(self, mi, m):
        K = _K(m)
        u = np.array([0.25*K, 0.5*K, 0.75*K, K])
        _, _, Ju = elliptic.jacobiEDJ(u, m, np.full(4, 0.3))
        np.testing.assert_allclose(np.asarray(Ju), _EXPECTED_JU[mi], atol=1e-12,
                                   err_msg=f"Ju table row m={m} vs MATLAB")


# ---------------------------------------------------------------------------
# Carlson — RC / RF / RD / RJ hardcoded tables
# Taken from matlab/tests/testCarlson.m (Group 5)
# ---------------------------------------------------------------------------

class TestCarlsonTable:
    def test_RC_exact_analytical(self):
        """DLMF analytical special values: RC(0,1)=pi/2, RC(0,1/4)=pi, etc."""
        x = np.array([0.0,  0.0,   0.25, 0.5,            1.0, 2.25])
        y = np.array([1.0,  0.25,  0.25, 1.0,            1.0, 2.0 ])
        expected = np.array([np.pi/2, np.pi, 2.0, np.pi/(2*np.sqrt(2)), 1.0, np.log(2)])
        RC = np.array([float(elliptic.carlsonRC(xi, yi)) for xi, yi in zip(x, y)])
        np.testing.assert_allclose(RC, expected, atol=1e-13, err_msg="RC analytical table vs MATLAB")

    def test_RC_generic(self):
        x = np.array([0.1,               1.5,               3.0              ])
        y = np.array([0.5,               0.3,               0.5              ])
        expected = np.array([1.750555828382159, 1.317852857616873, 0.9768180523022534])
        RC = np.array([float(elliptic.carlsonRC(xi, yi)) for xi, yi in zip(x, y)])
        np.testing.assert_allclose(RC, expected, atol=1e-13, err_msg="RC generic table vs MATLAB")

    def test_RF_table(self):
        xs = [1,   0.5, 0,   2,   0.25]
        ys = [2,   1,   0.5, 3,   0.5 ]
        zs = [3,   1.5, 1,   4,   1   ]
        expected = [0.7269459354689083, 1.028056801052127, 1.854074677301372,
                    0.5840828416771517, 1.370171633266872]
        for x, y, z, ref in zip(xs, ys, zs, expected):
            val = float(elliptic.carlsonRF(x, y, z))
            assert abs(val - ref) < 1e-12, f"RF({x},{y},{z})={val:.16g} expected {ref:.16g}"

    def test_RF_equals_K(self):
        """RF(0, 1-m, 1) = K(m) for m in {0.1..0.9} — table from testCarlson.m"""
        from scipy.special import ellipk
        m = np.array([0.1, 0.3, 0.5, 0.7, 0.9])
        expected_K = ellipk(m)
        for mi, ref in zip(m, expected_K):
            val = float(elliptic.carlsonRF(0.0, float(1 - mi), 1.0))
            assert abs(val - ref) < 1e-12, f"RF(0,1-{mi},1) vs K({mi})"

    def test_RD_table(self):
        xs = [0,   0,   0,   1,   0.5]
        ys = [2,   1,   0.5, 2,   1  ]
        zs = [1,   2,   1,   3,   1.5]
        expected = [1.797210352103389, 1.067937989667396, 3.020584777522179,
                    0.2904602810289907, 0.8215457375237986]
        for x, y, z, ref in zip(xs, ys, zs, expected):
            val = float(elliptic.carlsonRD(x, y, z))
            assert abs(val - ref) < 1e-12, f"RD({x},{y},{z})={val:.16g} expected {ref:.16g}"

    def test_RD_equals_D(self):
        """RD(0, 1-m, 1)/3 = D(m) = (K-E)/m — table from testCarlson.m"""
        from scipy.special import ellipk, ellipe
        m = np.array([0.1, 0.3, 0.5, 0.7, 0.9])
        expected_D = (ellipk(m) - ellipe(m)) / m
        for mi, ref in zip(m, expected_D):
            val = float(elliptic.carlsonRD(0.0, float(1 - mi), 1.0)) / 3.0
            assert abs(val - ref) < 1e-12, f"RD(0,1-{mi},1)/3 vs D({mi})"

    def test_RJ_table(self):
        xs = [0,   0,   1,   0.5]
        ys = [1,   1,   2,   1  ]
        zs = [2,   2,   3,   1.5]
        ps = [3,   0.5, 4,   2  ]
        expected = [0.7768862377858351, 2.936671269238118,
                    0.2398480997495678, 0.6783928711505076]
        for x, y, z, p, ref in zip(xs, ys, zs, ps, expected):
            val = float(elliptic.carlsonRJ(x, y, z, p))
            assert abs(val - ref) < 1e-12, f"RJ({x},{y},{z},{p})={val:.16g} expected {ref:.16g}"

    def test_RJ_equals_RD(self):
        """RJ(x,y,z,z) = RD(x,y,z) — DLMF identity"""
        for x, y, z in [(1.0, 2.0, 3.0), (0.5, 1.0, 1.5)]:
            rj = float(elliptic.carlsonRJ(x, y, z, z))
            rd = float(elliptic.carlsonRD(x, y, z))
            assert abs(rj - rd) < 1e-10, f"RJ({x},{y},{z},{z}) != RD({x},{y},{z}): {rj} vs {rd}"


# ---------------------------------------------------------------------------
# Bulirsch cel — generic reference table
# Taken from matlab/tests/testBulirschCEL.m (Group 4)
# ---------------------------------------------------------------------------

_CEL_TESTS = [
    (0.5, 0.3, 1.0, 0.5,  2.785372539236548),
    (0.5, 1.0, 0.8, 1.2,  2.229457648629679),
    (0.5, 2.0, 0.5, 1.5,  1.436498488170953),
    (0.7, 0.4, 1.0, 1.0,  3.066196685396346),
    (0.7, 1.0, 1.0, 1.0,  1.862640802332739),
    (0.7, 1.5, 0.3, 0.7,  0.7431229068315112),
    (0.3, 0.2, 0.5, 1.0,  6.654375216945862),
    (0.3, 1.0, 1.0, 0.09, 1.096477517392227),
    (0.9, 0.1, 1.0, 1.0,  5.378472453168364),
    (0.9, 1.0, 1.0, 1.0,  1.654616667522527),
]


class TestBulirschTable:
    @pytest.mark.parametrize("kc,p,a,b,ref", _CEL_TESTS,
                             ids=[f"cel(kc={t[0]},p={t[1]},a={t[2]},b={t[3]})" for t in _CEL_TESTS])
    def test_cel_generic(self, kc, p, a, b, ref):
        val = float(elliptic.cel(kc, p, a, b))
        assert abs(val - ref) < 1e-12, f"cel({kc},{p},{a},{b})={val:.16g} expected {ref:.16g}"

    def test_cel2_D(self):
        """cel2(kc, 0, 1) = D(m) — from testBulirschCEL.m"""
        m = np.array([0.1, 0.3, 0.5, 0.7, 0.9])
        kc = np.sqrt(1 - m)
        _, D, _ = elliptic.ellipticBD(m)
        C2D = np.array([float(elliptic.cel2(float(k), 0.0, 1.0)) for k in kc])
        np.testing.assert_allclose(C2D, np.asarray(D), atol=1e-12, err_msg="cel2(kc,0,1) != D(m)")

    def test_cel3_vs_Pi(self):
        """cel3(kc, 1-n) = Pi(n|m) — from testBulirschCEL.m"""
        m, n = 0.5, 0.3
        kc = float(np.sqrt(1 - m))
        Pi_leg = float(elliptic.elliptic3(np.pi / 2, m, n))
        Pi_cel = float(elliptic.cel3(kc, 1 - n))
        assert abs(Pi_cel - Pi_leg) < 1e-12, f"cel3 != Pi: {Pi_cel} vs {Pi_leg}"

    def test_cel_linearity(self):
        """cel(kc,p,a1+a2,b1+b2) = cel(kc,p,a1,b1) + cel(kc,p,a2,b2)"""
        kc, p = 0.6, 0.7
        c1 = float(elliptic.cel(kc, p, 1.0, 2.0))
        c2 = float(elliptic.cel(kc, p, 3.0, 0.5))
        c3 = float(elliptic.cel(kc, p, 4.0, 2.5))
        assert abs(c1 + c2 - c3) < 1e-12, f"cel not linear in a,b: {c1}+{c2}={c1+c2} != {c3}"
