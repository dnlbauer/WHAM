pub fn mean(x: &[f64]) -> f64 {
    x.iter().sum::<f64>() / x.len() as f64
}

pub fn autocov(x: &[f64]) -> f64 {
    let x_mean = mean(x);

    x.iter().map(|xi| {
        (xi-x_mean).powi(2)
    }).sum::<f64>() / x.len() as f64
}

pub fn sd(x: &[f64]) -> f64 {
    let x_mean = mean(x);
    let n = x.len() as f64;

    let sum = x.iter().map(|xi| {
        (xi-x_mean).powi(2)
    }).sum::<f64>();
    (1.0/(n-1.0) * sum).sqrt()
}

#[cfg(test)]
mod tests {
    use assert_approx_eq::assert_approx_eq;

    // a sine wave
    fn dataset() -> Vec<f64> {
        (0..100).map(|i| (i as f64 / 100.0 * std::f64::consts::PI).sin())
            .collect::<Vec<f64>>()
    }

    #[test]
    fn mean() {
        let ds = dataset();
        let expected =0.6366;
        
        let m = super::mean(&ds);
        assert_approx_eq!(m, expected, 0.0001);
    }

    #[test]
    fn autocorr() {
        let ds = dataset();
        let expected = 0.094_782;
        
        let m = super::autocov(&ds);
        assert_approx_eq!(m, expected, 0.000_001);
    }

    #[test]
    fn sd() {
        let ds = dataset();
        let expected = 0.309_418;
        
        let m = super::sd(&ds);
        assert_approx_eq!(m, expected, 0.000_001);
    }
}