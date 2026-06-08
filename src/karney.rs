//! Exact Transverse Mercator projection using Karney's Krüger n-series (order n^6).
//!
//! This is an alternative to the truncated OS Redfearn series used elsewhere in the
//! crate. The Redfearn series is expanded in powers of the longitude/easting offset
//! from the central meridian and is truncated at l^6 / e^7; its truncation error grows
//! steeply far from the -2 degree central meridian (a few mm at the extreme western
//! isles such as St Kilda). The Krüger n-series used here is expanded in the third
//! flattening `n` and is accurate to a few nanometres across the whole British National
//! Grid extent.
//!
//! Reference: C. F. F. Karney (2011), "Transverse Mercator with an accuracy of a few
//! nanometers", Journal of Geodesy 85(8), 475-485. The forward/inverse series and
//! conformal-latitude helpers follow the formulation used in GeographicLib.

/// tan of the conformal latitude from tan of the geographic latitude.
///
/// `e` is the (first) eccentricity of the ellipsoid.
fn taupf(tau: f64, e: f64) -> f64 {
    let tau1 = tau.hypot(1.0);
    let sig = (e * (e * tau / tau1).atanh()).sinh();
    sig.hypot(1.0) * tau - sig * tau1
}

/// tan of the geographic latitude from tan of the conformal latitude (inverse of [`taupf`]).
///
/// Newton's method; converges to full f64 precision in a handful of iterations because
/// the eccentricity of terrestrial ellipsoids is small.
fn tauf(taup: f64, e: f64) -> f64 {
    let e2m = 1.0 - e * e;
    // Initial guess.
    let mut tau = taup / e2m;
    for _ in 0..5 {
        let tau1 = tau.hypot(1.0);
        let sig = (e * (e * tau / tau1).atanh()).sinh();
        let taupa = sig.hypot(1.0) * tau - sig * tau1;
        let dtau = (taup - taupa) * (1.0 + e2m * tau * tau) / (e2m * tau1 * taupa.hypot(1.0));
        tau += dtau;
        if dtau.abs() < 1e-16 {
            break;
        }
    }
    tau
}

/// A Transverse Mercator projection parameterised by ellipsoid and grid constants,
/// evaluated with the Krüger n-series to order n^6.
pub(crate) struct TransverseMercator {
    a1f0: f64,     // rectifying radius * central scale factor
    e: f64,        // first eccentricity
    lon0: f64,     // central meridian (radians)
    fe: f64,       // false easting
    fn_: f64,      // false northing
    y0: f64,       // projected meridian distance to the true-origin latitude
    alp: [f64; 7], // forward series coefficients (1..=6 used)
    bet: [f64; 7], // inverse series coefficients (1..=6 used)
}

impl TransverseMercator {
    /// Build a projection from ellipsoid axes and grid parameters.
    ///
    /// * `a`, `b` - ellipsoid semi-major / semi-minor axes
    /// * `f0` - central meridian scale factor
    /// * `lat0`, `lon0` - true origin latitude / central meridian (radians)
    /// * `fe`, `fn_` - false easting / false northing
    pub(crate) fn new(a: f64, b: f64, f0: f64, lat0: f64, lon0: f64, fe: f64, fn_: f64) -> Self {
        let n = (a - b) / (a + b);
        let n2 = n * n;

        // Forward series coefficients (Karney 2011, eq. 35), Horner form in n.
        let alp = [
            0.0,
            n * (1. / 2.
                + n * (-2. / 3.
                    + n * (5. / 16.
                        + n * (41. / 180. + n * (-127. / 288. + n * (7891. / 37800.)))))),
            n2 * (13. / 48.
                + n * (-3. / 5.
                    + n * (557. / 1440. + n * (281. / 630. + n * (-1983433. / 1935360.))))),
            n2 * n
                * (61. / 240.
                    + n * (-103. / 140. + n * (15061. / 26880. + n * (167603. / 181440.)))),
            n2 * n2 * (49561. / 161280. + n * (-179. / 168. + n * (6601661. / 7257600.))),
            n2 * n2 * n * (34729. / 80640. + n * (-3418889. / 1995840.)),
            n2 * n2 * n2 * (212378941. / 319334400.),
        ];

        // Inverse series coefficients (Karney 2011, eq. 36), Horner form in n.
        let bet = [
            0.0,
            n * (1. / 2.
                + n * (-2. / 3.
                    + n * (37. / 96.
                        + n * (-1. / 360. + n * (-81. / 512. + n * (96199. / 604800.)))))),
            n2 * (1. / 48.
                + n * (1. / 15.
                    + n * (-437. / 1440. + n * (46. / 105. + n * (-1118711. / 3870720.))))),
            n2 * n * (17. / 480. + n * (-37. / 840. + n * (-209. / 4480. + n * (5569. / 90720.)))),
            n2 * n2 * (4397. / 161280. + n * (-11. / 504. + n * (-830251. / 7257600.))),
            n2 * n2 * n * (4583. / 161280. + n * (-108847. / 3991680.)),
            n2 * n2 * n2 * (20648693. / 638668800.),
        ];

        // Rectifying radius (Karney 2011, eq. 14), scaled by the central meridian factor.
        let a1 = a / (1. + n) * (1. + n2 * (1. / 4. + n2 * (1. / 64. + n2 / 256.)));
        let a1f0 = a1 * f0;

        let e2 = (a * a - b * b) / (a * a);
        let e = e2.sqrt();

        let mut proj = TransverseMercator {
            a1f0,
            e,
            lon0,
            fe,
            fn_,
            y0: 0.0,
            alp,
            bet,
        };
        // Projected meridian distance from the equator to the true-origin latitude, so
        // that the true origin (lat0, lon0) maps to (fe, fn_).
        proj.y0 = proj.meridian(lat0);
        proj
    }

    /// Projected northing (relative to the equator) of a point on the central meridian.
    fn meridian(&self, lat: f64) -> f64 {
        let taup = taupf(lat.tan(), self.e);
        let xip = taup.atan(); // atan2(taup, 1) since eta' = 0 on the central meridian
        let mut xi = xip;
        for p in 1..=6 {
            xi += self.alp[p] * ((2 * p) as f64 * xip).sin();
        }
        self.a1f0 * xi
    }

    /// Forward projection: (longitude, latitude) in radians to (easting, northing).
    pub(crate) fn forward(&self, lon: f64, lat: f64) -> (f64, f64) {
        let lam = lon - self.lon0;
        let taup = taupf(lat.tan(), self.e);
        let (sin_lam, cos_lam) = lam.sin_cos();
        let xip = taup.atan2(cos_lam);
        let etap = (sin_lam / taup.hypot(cos_lam)).asinh();

        let mut xi = xip;
        let mut eta = etap;
        for p in 1..=6 {
            let s = (2 * p) as f64;
            let (sin_s_xip, cos_s_xip) = (s * xip).sin_cos();
            xi += self.alp[p] * sin_s_xip * (s * etap).cosh();
            eta += self.alp[p] * cos_s_xip * (s * etap).sinh();
        }

        let easting = self.fe + self.a1f0 * eta;
        let northing = self.fn_ + (self.a1f0 * xi - self.y0);
        (easting, northing)
    }

    /// Inverse projection: (easting, northing) to (longitude, latitude) in radians.
    pub(crate) fn inverse(&self, easting: f64, northing: f64) -> (f64, f64) {
        let xi = (northing - self.fn_ + self.y0) / self.a1f0;
        let eta = (easting - self.fe) / self.a1f0;

        let mut xip = xi;
        let mut etap = eta;
        for p in 1..=6 {
            let s = (2 * p) as f64;
            let (sin_s_xi, cos_s_xi) = (s * xi).sin_cos();
            xip -= self.bet[p] * sin_s_xi * (s * eta).cosh();
            etap -= self.bet[p] * cos_s_xi * (s * eta).sinh();
        }

        let sinh_etap = etap.sinh();
        let cos_xip = xip.cos();
        let taup = xip.sin() / sinh_etap.hypot(cos_xip);
        let lam = sinh_etap.atan2(cos_xip);
        let lon = self.lon0 + lam;
        let lat = tauf(taup, self.e).atan();
        (lon, lat)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const A: f64 = 6378137.000;
    const B: f64 = 6356752.3141;
    const F0: f64 = 0.9996012717;
    const PI: f64 = std::f64::consts::PI;

    fn bng() -> TransverseMercator {
        TransverseMercator::new(
            A,
            B,
            F0,
            49.0 * PI / 180.,
            -2.0 * PI / 180.,
            400000.,
            -100000.,
        )
    }

    #[test]
    fn taup_round_trip() {
        let e = ((A * A - B * B) / (A * A)).sqrt();
        for deg in [-60.0, -10.0, 0.0, 1.0, 49.0, 52.658, 60.84] {
            let tau = (deg * PI / 180.).tan();
            let back = tauf(taupf(tau, e), e);
            assert!((tau - back).abs() <= 1e-12 * (1.0 + tau.abs()), "deg={deg}");
        }
    }

    #[test]
    fn caister_reference() {
        // OSGM15 user guide p20-23: lon 1.716073973, lat 52.658007833 -> ETRS89 (651307.003, 313255.686)
        let p = bng();
        let (e, n) = p.forward(1.716073973_f64.to_radians(), 52.658007833_f64.to_radians());
        assert!((e - 651307.003).abs() < 1e-3, "easting {e}");
        assert!((n - 313255.686).abs() < 1e-3, "northing {n}");
    }

    #[test]
    fn forward_inverse_round_trip() {
        let p = bng();
        // St Kilda area (far west, where the Redfearn series struggles)
        let cases: [(f64, f64); 3] = [
            (-8.578544561, 57.81351838),
            (-2.0, 49.0),
            (1.716073973, 52.658007833),
        ];
        for (lon, lat) in cases {
            let (e, n) = p.forward(lon.to_radians(), lat.to_radians());
            let (lon2, lat2) = p.inverse(e, n);
            assert!(
                (lon - lon2.to_degrees()).abs() < 1e-9,
                "lon {lon} -> {}",
                lon2.to_degrees()
            );
            assert!(
                (lat - lat2.to_degrees()).abs() < 1e-9,
                "lat {lat} -> {}",
                lat2.to_degrees()
            );
        }
    }
}
