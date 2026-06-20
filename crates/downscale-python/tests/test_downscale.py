"""Smoke tests de los bindings Python de downscale-rs.

Ejecutar desde la raíz del workspace con el wheel instalado:
    .venv/bin/python crates/downscale-python/tests/test_downscale.py
"""

import numpy as np

import downscale_rs as ds

rng = np.random.default_rng(42)


def test_eqm_removes_additive_bias():
    obs = rng.normal(10.0, 2.0, 2000)
    model = obs + 3.0
    qm = ds.QuantileMapping(obs, model, n_quantiles=100, kind="add")
    corrected = qm.apply(model)
    assert abs(corrected.mean() - obs.mean()) < 1e-6
    assert abs(qm.correct_one(13.0) - 10.0) < 0.2


def test_eqm_midpoint_nodes_and_errors():
    obs = rng.normal(0.0, 1.0, 500)
    ds.QuantileMapping(obs, obs, nodes="midpoint")
    try:
        ds.QuantileMapping(obs, obs, kind="bogus")
        raise AssertionError("kind inválido no rechazado")
    except ValueError:
        pass
    try:
        ds.QuantileMapping(np.array([1.0, np.nan]), obs)
        raise AssertionError("NaN no rechazado")
    except ValueError as e:
        assert "no finito" in str(e)


def test_parametric_gamma_corrects_wet_day_frequency():
    u = rng.uniform(size=6000)
    obs = np.where(u < 0.7, 0.0, rng.gamma(0.8, 8.0, 6000))
    v = rng.uniform(size=6000)
    model = np.where(v < 0.4, 0.0, rng.gamma(0.8, 5.0, 6000))
    pqm = ds.ParametricQuantileMapping(obs, model, dist="gamma", wet_threshold=0.1)
    corrected = pqm.apply(model)
    dry = lambda s: float(np.mean(s < 0.1))
    assert abs(dry(corrected) - dry(obs)) < 0.03
    assert corrected.min() >= 0.0


def test_delta_change():
    dc = ds.DeltaChange(np.array([2.0, 4.0]), np.array([3.0, 6.0]), kind="mult")
    assert dc.delta == 1.5
    np.testing.assert_allclose(dc.apply(np.array([10.0, 0.0])), [15.0, 0.0])


def test_wet_day_correction():
    obs = np.concatenate([np.zeros(80), np.arange(1.0, 21.0)])
    model = np.linspace(0.1, 40.0, 100)
    wd = ds.WetDayCorrection(obs, model, obs_wet_threshold=0.1)
    transformed = wd.transform(model)
    assert 78 <= int(np.sum(transformed == 0.0)) <= 82
    assert wd.model_threshold > 0.0


def test_analog_and_regression():
    x = np.linspace(0.0, 5.0, 500).reshape(-1, 1)
    y = np.sin(x[:, 0])
    ad = ds.AnalogDownscaling(x, y, k=3)
    pred = ad.predict(np.array([[1.7], [3.3]]))
    np.testing.assert_allclose(pred, np.sin([1.7, 3.3]), atol=0.02)

    x2 = rng.uniform(size=(800, 2)) * [10.0, 4.0]
    y2 = 2.0 + 3.0 * x2[:, 0] - x2[:, 1]
    lm = ds.LinearDownscaling(x2, y2)
    assert abs(lm.intercept - 2.0) < 1e-9
    np.testing.assert_allclose(lm.coefs, [3.0, -1.0], atol=1e-9)
    assert lm.r2 > 0.999999


def test_metrics():
    a = np.array([0.0, 0.0])
    b = np.array([3.0, 4.0])
    assert abs(ds.rmse(a, b) - np.sqrt(12.5)) < 1e-12
    assert ds.mean_bias(b, a) == 3.5
    assert ds.ks_statistic(a, a) == 0.0
    qb = ds.quantile_bias(b, a, probs=np.array([0.5]))
    assert qb[0][0] == 0.5 and qb[0][3] == 3.5


def test_validate_split_dict():
    obs = 15.0 + 8.0 * np.sin(np.arange(730) * 2 * np.pi / 365)
    model = obs * 1.1 + 3.0
    report = ds.validate_split(obs, model, calib_frac=0.5, kind="add")
    assert report["split_index"] == 365
    assert report["rmse"] < report["rmse_raw"]
    assert report["ks"] <= report["ks_raw"]
    assert len(report["quantile_bias"]) == 12


def test_qdm_preserves_change_signal():
    obs = 10.0 + 8.0 * rng.uniform(size=2000)
    hist = obs + 2.5  # sesgo
    proj = hist + 4.0  # señal de cambio
    qdm = ds.QuantileDeltaMapping(obs, hist, n_quantiles=100, kind="add")
    corrected = qdm.apply(proj)
    assert abs(corrected.mean() - (obs.mean() + 4.0)) < 1e-6


def test_schaake_shuffle_restores_rank_structure():
    template = rng.uniform(size=(300, 2))
    template[:, 1] = template[:, 0] + 0.01 * rng.uniform(size=300)  # corr alta
    corrected = rng.uniform(size=(300, 2)) * 50.0
    out = ds.schaake_shuffle(template, corrected)
    # Marginales preservadas y correlación de rangos alta como la plantilla.
    np.testing.assert_allclose(np.sort(out[:, 0]), np.sort(corrected[:, 0]))
    from numpy import corrcoef

    rank = lambda x: np.argsort(np.argsort(x))
    rho = corrcoef(rank(out[:, 0]), rank(out[:, 1]))[0, 1]
    assert rho > 0.99, f"rho = {rho}"


def test_hargreaves_fao56():
    # Santiago en verano: PET plausible.
    pet = ds.hargreaves(np.array([13.0]), np.array([30.0]), [15], -33.45)
    assert 4.0 <= pet[0] <= 8.0


def test_qdm_windowed():
    obs = 10.0 + 5.0 * rng.uniform(size=2000)
    hist = obs + 2.0
    qdm = ds.QuantileDeltaMapping(obs, hist, n_quantiles=100, kind="add")
    # Proyección no estacionaria: dos regímenes de cambio.
    proj = np.concatenate([hist[:1000] + 1.0, hist[1000:] + 9.0])
    w = qdm.apply_windowed(proj, 1000)
    assert abs(w[:1000].mean() - (obs.mean() + 1.0)) < 0.2
    assert abs(w[1000:].mean() - (obs.mean() + 9.0)) < 0.2


def test_mbcn_recovers_dependence():
    n = 3000
    z1, z2 = rng.standard_normal(n), rng.standard_normal(n)
    w1, w2 = rng.standard_normal(n), rng.standard_normal(n)
    obs = np.column_stack([10 + 2 * z1, 20 + 0.8 * 3 * z1 + 0.6 * 3 * z2])
    model = np.column_stack([13 + 3 * w1, 17 - 0.3 * 4 * w1 + 0.95 * 4 * w2])
    out = ds.mbcn(obs, model, n_iterations=30, seed=42)
    corr = lambda m: np.corrcoef(m[:, 0], m[:, 1])[0, 1]
    # La correlación corregida se acerca a la observada.
    assert abs(corr(out) - corr(obs)) < 0.12
    # Determinismo.
    out2 = ds.mbcn(obs, model, seed=42)
    assert np.array_equal(out, out2)


def test_correct_grid_3d():
    # Campo [time, lat, lon] con sesgo multiplicativo por celda + máscara de mar.
    obs = rng.gamma(1.0, 3.0, (500, 2, 3))
    factor = 1.0 + 0.2 * np.arange(6).reshape(2, 3)
    model = obs * factor
    model[:, 0, 0] = np.nan  # mar
    out = ds.correct_grid(obs, model, kind="mult")
    assert out.shape == model.shape
    assert np.isnan(out[:, 0, 0]).all()  # celda enmascarada
    # celda válida: sesgo medio corregido ≈ 0
    assert abs(out[:, 1, 2].mean() - obs[:, 1, 2].mean()) < 1e-6


def test_real_data_matches_rust_cli():
    """Paridad bindings vs CLI sobre el caso Quinta Normal (si hay datos)."""
    import os

    base = os.path.join(os.path.dirname(__file__), "../../../data/parity")
    if not os.path.isdir(base):
        import pytest

        pytest.skip("data/parity no disponible")
    load = lambda n: np.loadtxt(
        os.path.join(base, f"{n}.csv"), delimiter=",", skiprows=1, usecols=1
    )
    obs_cal, mod_cal = load("obs_cal"), load("model_cal")
    mod_val, rust_out = load("model_val"), load("rust_corrected")
    qm = ds.QuantileMapping(obs_cal, mod_cal, n_quantiles=100, kind="mult")
    corrected = qm.apply(mod_val)
    # El CSV del CLI guarda 3 decimales.
    assert np.max(np.abs(corrected - rust_out)) < 6e-4


if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    passed = skipped = 0
    for t in tests:
        try:
            t()
            print(f"✓ {t.__name__}")
            passed += 1
        except Exception as e:  # pytest.skip lanza Skipped fuera de pytest
            if type(e).__name__ == "Skipped":
                print(f"– {t.__name__} (omitido: {e})")
                skipped += 1
            else:
                raise
    print(f"\n{passed} OK, {skipped} omitidos")
