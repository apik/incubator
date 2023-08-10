#![allow(non_snake_case)]

use crate::integrator::IntegrandFunction;
use crate::integrator::HyperCube;


use intbits::Bits;

// (a,0,0,0,0,...)
// with points saved for bisection direction calculation
// pub fn orbitAdiff(lambdaA:f64, hc:&HyperCube, f:IntegrandFunction, n:usize, lPlus:& mut Vec<(f64,f64)>, lMinus:& mut Vec<f64>) -> f64 {
pub fn orbitAdiff(lambda_a:f64, hc:&HyperCube, f:IntegrandFunction, n:usize, l_shifts: & mut [(f64,f64)], f_evals:&mut usize) -> f64 {
    let mut res: f64 = 0.0;
    // Position of the first generator
    for g_a_pos in 0..n {
        for sgns in 0..2_i32.pow(1) {
            let mut x = hc.xc.to_vec();
            if sgns.bit(0) {x[g_a_pos] += lambda_a*hc.hw[g_a_pos]} else {x[g_a_pos] -= lambda_a*hc.hw[g_a_pos]}
            let f_val = f(&x);
            *f_evals += 1;      // increase function evaluations counter
            // Store function values first - "+" "shift, second "-" shift
            if sgns.bit(0) {l_shifts[g_a_pos].0 = f_val} else {l_shifts[g_a_pos].1 = f_val};
            res += f_val;
        }
    }
    res*hc.hwvol
}

// (a,0,0,0,0,...)
pub fn orbitA(lambda_a:f64, hc:&HyperCube, f:IntegrandFunction, n:usize, f_evals:&mut usize) -> f64 {
    let mut res: f64 = 0.0;
    // Position of the first generator
    for g_a_pos in 0..n {
        for sgns in 0..2_i32.pow(1) {
            let mut x = hc.xc.to_vec();
            if sgns.bit(0) {x[g_a_pos] += lambda_a*hc.hw[g_a_pos]} else {x[g_a_pos] -= lambda_a*hc.hw[g_a_pos]}
            res += f(&x);
            *f_evals += 1;
        }
    }
    res*hc.hwvol
}

// (a,a,0,0,0,...)
pub fn orbitAA(lambda_a:f64, hc:&HyperCube, f:IntegrandFunction, n:usize, f_evals:&mut usize) -> f64 {
    let mut res: f64 = 0.0;
    // Position of the first generator of the A type
    for gA1pos in 0..n {
        for gA2pos in (gA1pos+1)..n {
            for sgns in 0..2_i32.pow(2) {
                let mut x = hc.xc.to_vec();
                if sgns.bit(0) {x[gA1pos] += lambda_a*hc.hw[gA1pos]} else {x[gA1pos] -= lambda_a*hc.hw[gA1pos]}
                if sgns.bit(1) {x[gA2pos] += lambda_a*hc.hw[gA2pos]} else {x[gA2pos] -= lambda_a*hc.hw[gA2pos]}
                res += f(&x);
                *f_evals += 1;
            }
        }}
    res*hc.hwvol
}

// (a,a,a,0,0,...)
pub fn orbitAAA(lambda_a:f64, hc:&HyperCube, f:IntegrandFunction, n:usize, f_evals:&mut usize) -> f64 {
    let mut res: f64 = 0.0;
    // Position of the first generator of the A type
    for gA1pos in 0..n {
        for gA2pos in (gA1pos+1)..n {
            for gA3pos in (gA2pos+1)..n {
                for sgns in 0..2_i32.pow(3) {
                    let mut x = hc.xc.to_vec();
                    if sgns.bit(0) {x[gA1pos] += lambda_a*hc.hw[gA1pos]} else {x[gA1pos] -= lambda_a*hc.hw[gA1pos]}
                    if sgns.bit(1) {x[gA2pos] += lambda_a*hc.hw[gA2pos]} else {x[gA2pos] -= lambda_a*hc.hw[gA2pos]}
                    if sgns.bit(2) {x[gA3pos] += lambda_a*hc.hw[gA3pos]} else {x[gA3pos] -= lambda_a*hc.hw[gA3pos]}
                    res += f(&x);
                    *f_evals += 1;
                }
            }}}
    res*hc.hwvol
}

// (a,a,a,a,0,...)
pub fn orbitAAAA(lambda_a:f64, hc:&HyperCube, f:IntegrandFunction, n:usize, f_evals:&mut usize) -> f64 {
    let mut res: f64 = 0.0;
    // Position of the first generator of the A type
    for gA1pos in 0..n {
        for gA2pos in (gA1pos+1)..n {
            for gA3pos in (gA2pos+1)..n {
                for gA4pos in (gA3pos+1)..n {
                    for sgns in 0..2_i32.pow(4) {
                        let mut x = hc.xc.to_vec();
                        if sgns.bit(0) {x[gA1pos] += lambda_a*hc.hw[gA1pos]} else {x[gA1pos] -= lambda_a*hc.hw[gA1pos]}
                        if sgns.bit(1) {x[gA2pos] += lambda_a*hc.hw[gA2pos]} else {x[gA2pos] -= lambda_a*hc.hw[gA2pos]}
                        if sgns.bit(2) {x[gA3pos] += lambda_a*hc.hw[gA3pos]} else {x[gA3pos] -= lambda_a*hc.hw[gA3pos]}
                        if sgns.bit(3) {x[gA4pos] += lambda_a*hc.hw[gA4pos]} else {x[gA4pos] -= lambda_a*hc.hw[gA4pos]}
                        res += f(&x);
                        *f_evals += 1;
                    }
                }}}}
    res*hc.hwvol
}

// (a,a,a,a,a,...)
pub fn orbitAAAAA(lambda_a:f64, hc:&HyperCube, f:IntegrandFunction, n:usize, f_evals:&mut usize) -> f64 {
    let mut res: f64 = 0.0;
    // Position of the first generator of the A type
    for gA1pos in 0..n {
        for gA2pos in (gA1pos+1)..n {
            for gA3pos in (gA2pos+1)..n {
                for gA4pos in (gA3pos+1)..n {
                    for gA5pos in (gA4pos+1)..n {
                        for sgns in 0..2_i32.pow(5) {
                            let mut x = hc.xc.to_vec();
                            if sgns.bit(0) {x[gA1pos] += lambda_a*hc.hw[gA1pos]} else {x[gA1pos] -= lambda_a*hc.hw[gA1pos]}
                            if sgns.bit(1) {x[gA2pos] += lambda_a*hc.hw[gA2pos]} else {x[gA2pos] -= lambda_a*hc.hw[gA2pos]}
                            if sgns.bit(2) {x[gA3pos] += lambda_a*hc.hw[gA3pos]} else {x[gA3pos] -= lambda_a*hc.hw[gA3pos]}
                            if sgns.bit(3) {x[gA4pos] += lambda_a*hc.hw[gA4pos]} else {x[gA4pos] -= lambda_a*hc.hw[gA4pos]}
                            if sgns.bit(4) {x[gA5pos] += lambda_a*hc.hw[gA5pos]} else {x[gA5pos] -= lambda_a*hc.hw[gA5pos]}
                            res += f(&x);
                            *f_evals += 1;
                        }
                    }}}}}
    res*hc.hwvol
}

// (b,a,0,0,0,...)
pub fn orbitBA(lambda_b:f64, lambda_a:f64, hc:&HyperCube, f:IntegrandFunction, n:usize, f_evals:&mut usize) -> f64 {
    let mut res: f64 = 0.0;
    // Loop over AB configuration
    for g_a_pos in 0..n {
        for gBpos in (g_a_pos+1)..n {
            for sgns in 0..2_i32.pow(2) {
                let mut x = hc.xc.to_vec();
                if sgns.bit(0) {x[g_a_pos] += lambda_a*hc.hw[g_a_pos]} else {x[g_a_pos] -= lambda_a*hc.hw[g_a_pos]}
                if sgns.bit(1) {x[gBpos] += lambda_b*hc.hw[gBpos]} else {x[gBpos] -= lambda_b*hc.hw[gBpos]}
                res += f(&x);
                *f_evals += 1;
            }
        }}
    // Loop over BA configuration
    for gBpos in 0..n {
        for g_a_pos in (gBpos+1)..n {
            for sgns in 0..2_i32.pow(2) {
                let mut x = hc.xc.to_vec();
                if sgns.bit(0) {x[g_a_pos] += lambda_a*hc.hw[g_a_pos]} else {x[g_a_pos] -= lambda_a*hc.hw[g_a_pos]}
                if sgns.bit(1) {x[gBpos] += lambda_b*hc.hw[gBpos]} else {x[gBpos] -= lambda_b*hc.hw[gBpos]}
                res += f(&x);
                *f_evals += 1;
            }
        }}
    res*hc.hwvol
}

// (b,a,a,0,0,...)
pub fn orbitBAA(lambda_b:f64, lambda_a:f64, hc:&HyperCube, f:IntegrandFunction, n:usize, f_evals:&mut usize) -> f64 {
    let mut res: f64 = 0.0;
    // Loop over AAB configuration
    for gA1pos in 0..n {
        for gA2pos in (gA1pos+1)..n {
            for gBpos in (gA2pos+1)..n {
                for sgns in 0..2_i32.pow(3) {
                    let mut x = hc.xc.to_vec();
                    if sgns.bit(0) {x[gA1pos] += lambda_a*hc.hw[gA1pos]} else {x[gA1pos] -= lambda_a*hc.hw[gA1pos]}
                    if sgns.bit(1) {x[gA2pos] += lambda_a*hc.hw[gA2pos]} else {x[gA2pos] -= lambda_a*hc.hw[gA2pos]}
                    if sgns.bit(2) {x[gBpos]  += lambda_b*hc.hw[gBpos]}  else {x[gBpos]  -= lambda_b*hc.hw[gBpos]}
                    res += f(&x);
                    *f_evals += 1;
                }
            }}}
    // Loop over ABA configuration
    for gA1pos in 0..n {
        for gBpos in (gA1pos+1)..n {
            for gA2pos in (gBpos+1)..n {
                for sgns in 0..2_i32.pow(3) {
                    let mut x = hc.xc.to_vec();
                    if sgns.bit(0) {x[gA1pos] += lambda_a*hc.hw[gA1pos]} else {x[gA1pos] -= lambda_a*hc.hw[gA1pos]}
                    if sgns.bit(1) {x[gA2pos] += lambda_a*hc.hw[gA2pos]} else {x[gA2pos] -= lambda_a*hc.hw[gA2pos]}
                    if sgns.bit(2) {x[gBpos]  += lambda_b*hc.hw[gBpos]}  else {x[gBpos]  -= lambda_b*hc.hw[gBpos]}
                    res += f(&x);
                    *f_evals += 1;
                }
            }}}
    // Loop over BAA configuration
    for gBpos in 0..n {
        for gA1pos in (gBpos+1)..n {
            for gA2pos in (gA1pos+1)..n {
                for sgns in 0..2_i32.pow(3) {
                    let mut x = hc.xc.to_vec();
                    if sgns.bit(0) {x[gA1pos] += lambda_a*hc.hw[gA1pos]} else {x[gA1pos] -= lambda_a*hc.hw[gA1pos]}
                    if sgns.bit(1) {x[gA2pos] += lambda_a*hc.hw[gA2pos]} else {x[gA2pos] -= lambda_a*hc.hw[gA2pos]}
                    if sgns.bit(2) {x[gBpos]  += lambda_b*hc.hw[gBpos]}  else {x[gBpos]  -= lambda_b*hc.hw[gBpos]}
                    res += f(&x);
                    *f_evals += 1;
                }
            }}}
    res*hc.hwvol
}

// (b,a,a,a,0,...)
pub fn orbitBAAA(lambda_b:f64, lambda_a:f64, hc:&HyperCube, f:IntegrandFunction, n:usize, f_evals:&mut usize) -> f64 {
    let mut res: f64 = 0.0;
    // Loop over AAAB configuration
    for gA1pos in 0..n {
        for gA2pos in (gA1pos+1)..n {
            for gA3pos in (gA2pos+1)..n {
                for gBpos in (gA3pos+1)..n {
                    for sgns in 0..2_i32.pow(4) {
                        let mut x = hc.xc.to_vec();
                        if sgns.bit(0) {x[gA1pos] += lambda_a*hc.hw[gA1pos]} else {x[gA1pos] -= lambda_a*hc.hw[gA1pos]}
                        if sgns.bit(1) {x[gA2pos] += lambda_a*hc.hw[gA2pos]} else {x[gA2pos] -= lambda_a*hc.hw[gA2pos]}
                        if sgns.bit(2) {x[gA3pos] += lambda_a*hc.hw[gA3pos]} else {x[gA3pos] -= lambda_a*hc.hw[gA3pos]}
                        if sgns.bit(3) {x[gBpos]  += lambda_b*hc.hw[gBpos]}  else {x[gBpos]  -= lambda_b*hc.hw[gBpos]}
                        res += f(&x);
                        *f_evals += 1;
                    }
                }}}}
    // Loop over AABA configuration
    for gA1pos in 0..n {
        for gA2pos in (gA1pos+1)..n {
            for gBpos in (gA2pos+1)..n {
                for gA3pos in (gBpos+1)..n {
                    for sgns in 0..2_i32.pow(4) {
                        let mut x = hc.xc.to_vec();
                        if sgns.bit(0) {x[gA1pos] += lambda_a*hc.hw[gA1pos]} else {x[gA1pos] -= lambda_a*hc.hw[gA1pos]}
                        if sgns.bit(1) {x[gA2pos] += lambda_a*hc.hw[gA2pos]} else {x[gA2pos] -= lambda_a*hc.hw[gA2pos]}
                        if sgns.bit(2) {x[gA3pos] += lambda_a*hc.hw[gA3pos]} else {x[gA3pos] -= lambda_a*hc.hw[gA3pos]}
                        if sgns.bit(3) {x[gBpos]  += lambda_b*hc.hw[gBpos]}  else {x[gBpos]  -= lambda_b*hc.hw[gBpos]}
                        res += f(&x);
                        *f_evals += 1;
                    }
                }}}}
    // Loop over ABAA configuration
    for gA1pos in 0..n {
        for gBpos in (gA1pos+1)..n {
            for gA2pos in (gBpos+1)..n {
                for gA3pos in (gA2pos+1)..n {
                    for sgns in 0..2_i32.pow(4) {
                        let mut x = hc.xc.to_vec();
                        if sgns.bit(0) {x[gA1pos] += lambda_a*hc.hw[gA1pos]} else {x[gA1pos] -= lambda_a*hc.hw[gA1pos]}
                        if sgns.bit(1) {x[gA2pos] += lambda_a*hc.hw[gA2pos]} else {x[gA2pos] -= lambda_a*hc.hw[gA2pos]}
                        if sgns.bit(2) {x[gA3pos] += lambda_a*hc.hw[gA3pos]} else {x[gA3pos] -= lambda_a*hc.hw[gA3pos]}
                        if sgns.bit(3) {x[gBpos]  += lambda_b*hc.hw[gBpos]}  else {x[gBpos]  -= lambda_b*hc.hw[gBpos]}
                        res += f(&x);
                        *f_evals += 1;
                    }
                }}}}
    // Loop over BAAA configuration
    for gBpos in 0..n {
        for gA1pos in (gBpos+1)..n {
            for gA2pos in (gA1pos+1)..n {
                for gA3pos in (gA2pos+1)..n {
                    for sgns in 0..2_i32.pow(4) {
                        let mut x = hc.xc.to_vec();
                        if sgns.bit(0) {x[gA1pos] += lambda_a*hc.hw[gA1pos]} else {x[gA1pos] -= lambda_a*hc.hw[gA1pos]}
                        if sgns.bit(1) {x[gA2pos] += lambda_a*hc.hw[gA2pos]} else {x[gA2pos] -= lambda_a*hc.hw[gA2pos]}
                        if sgns.bit(2) {x[gA3pos] += lambda_a*hc.hw[gA3pos]} else {x[gA3pos] -= lambda_a*hc.hw[gA3pos]}
                        if sgns.bit(3) {x[gBpos]  += lambda_b*hc.hw[gBpos]}  else {x[gBpos]  -= lambda_b*hc.hw[gBpos]}
                        res += f(&x);
                        *f_evals += 1;
                    }
                }}}}
    res*hc.hwvol
}

// (delta,delta,delta,delta,delta,...,delta)
pub fn orbitDelta(delta:f64, hc:&HyperCube, f:IntegrandFunction, n:usize, f_evals:&mut usize) -> f64 {
    let mut res: f64 = 0.0;
    // Position of the first generator
    for sgns in 0..2_i32.pow(n as u32) {
        let mut x = hc.xc.to_vec();
        for (s_pos, x_i) in x.iter_mut().enumerate() {
            if sgns.bit(s_pos) {*x_i += delta*hc.hw[s_pos]} else {*x_i -= delta*hc.hw[s_pos]}
        }
        // for s_pos in 0..n {
        //     if sgns.bit(s_pos) {x[s_pos] += delta*hc.hw[s_pos]} else {x[s_pos] -= delta*hc.hw[s_pos]}
        // }
        res += f(&x);
        *f_evals += 1;
    }
    res*hc.hwvol
}

#[cfg(test)]
mod tests {

    use crate::rules::Generators;
    use crate::integrator::Point;
    use assert_approx_eq::assert_approx_eq;
    use super::*;

    // Test functions
    fn f_11(x:&Point) -> f64 { x[0]*x[1] }
    fn f_12(x:&Point) -> f64 { x[0]*x[1].powi(2) }
    fn f_22(x:&Point) -> f64 { x[0].powi(2)*x[1].powi(2) }
    fn f_222(x:&Point) -> f64 { x[0].powi(2)*x[1].powi(2)*x[2].powi(2) }
    fn f_2222(x:&Point) -> f64 { x[0].powi(2)*x[1].powi(2)*x[2].powi(2)*x[3].powi(2) }
    fn f_1234(x:&Point) -> f64 { x[0].powi(1)*x[1].powi(2)*x[2].powi(3)*x[3].powi(4) }

    #[test]
    fn orbit_A() {
        let mut f_evals:usize = 0;
        // Initialize generators
        let gen = Generators::new(100);
        let cube5d = HyperCube::new(&vec![1.0; 5], &vec![1.0; 5]);
        const ERR:f64 = 1e-9f64;

        // (1,0,0,0,0,...)
        assert_approx_eq!(orbitA(gen.get_l1(), &cube5d, f_11 , 5, &mut f_evals), 10.0, ERR);
        assert_approx_eq!(orbitA(gen.get_l1(), &cube5d, f_12 , 5, &mut f_evals), 11.827517592817278, ERR);
        assert_approx_eq!(orbitA(gen.get_l1(), &cube5d, f_22 , 5, &mut f_evals), 13.655035185634554720, ERR);
        assert_approx_eq!(orbitA(gen.get_l1(), &cube5d, f_222 , 5, &mut f_evals), 15.482552778451832080, ERR);
        assert_approx_eq!(orbitA(gen.get_l1(), &cube5d, f_2222 , 5, &mut f_evals), 17.310070371269109439, ERR);

        // (2,0,0,0,0,...)
        assert_approx_eq!(orbitA(gen.get_l2(), &cube5d, f_2222 , 5, &mut f_evals), 11.319059433251210865, ERR);

        // (3,0,0,0,0,...)
        assert_approx_eq!(orbitA(gen.get_l3(), &cube5d, f_2222 , 5, &mut f_evals), 16.411847955508167496, ERR);

        // (4,0,0,0,0,...)
        assert_approx_eq!(orbitA(gen.get_l4(), &cube5d, f_2222 , 5, &mut f_evals), 15.745865464580127547, ERR);
    }

    #[test]
    fn orbit_AA() {
        let mut f_evals:usize = 0;
        // Initialize generators
        let gen = Generators::new(100);
        let cube5d = HyperCube::new(&vec![1.0; 5], &vec![1.0; 5]);
        const ERR:f64 = 1e-9f64;

        // (1,1,0,0,0,...)
        assert_approx_eq!(orbitAA(gen.get_l1(), &cube5d, f_1234 , 5, &mut f_evals), 301.94300617061102444, ERR);

        // (2,2,0,0,0,...)
        assert_approx_eq!(orbitAA(gen.get_l2(), &cube5d, f_1234 , 5, &mut f_evals), 69.823999848017142088, ERR);

    }

    #[test]
    fn orbit_AAA() {
        let mut f_evals:usize = 0;
        // Initialize generators
        let gen = Generators::new(100);
        let cube5d = HyperCube::new(&vec![1.0; 5], &vec![1.0; 5]);
        const ERR:f64 = 1e-9f64;

        // (1,1,1,0,0,...)
        assert_approx_eq!(orbitAAA(gen.get_l1(), &cube5d, f_1234 , 5, &mut f_evals), 1299.5720247860022532, ERR);

    }

    #[test]
    fn orbit_AAAA() {
        let mut f_evals:usize = 0;
        // Initialize generators
        let gen = Generators::new(100);
        let cube5d = HyperCube::new(&vec![1.0; 5], &vec![1.0; 5]);
        const ERR:f64 = 1e-9f64;

        // (1,1,1,1,0,...)
        assert_approx_eq!(orbitAAAA(gen.get_l1(), &cube5d, f_1234 , 5, &mut f_evals), 2453.2147836311119702, ERR);

    }

    #[test]
    fn orbit_AAAAA() {
        let mut f_evals:usize = 0;
        // Initialize generators
        let gen = Generators::new(100);
        let cube5d = HyperCube::new(&vec![1.0; 5], &vec![1.0; 5]);
        const ERR:f64 = 1e-9f64;

        // (1,1,1,1,1,...)
        assert_approx_eq!(orbitAAAAA(gen.get_l1(), &cube5d, f_1234 , 5, &mut f_evals), 1676.5641382158854978, ERR);

    }

    #[test]
    fn orbit_BA() {
        let mut f_evals:usize = 0;
        // Initialize generators
        let gen = Generators::new(100);
        let cube5d = HyperCube::new(&vec![1.0; 5], &vec![1.0; 5]);
        const ERR:f64 = 1e-9f64;

        // (2,1,0,0,0,...)
        assert_approx_eq!(orbitBA(gen.get_l2(), gen.get_l1(), &cube5d, f_1234 , 5, &mut f_evals), 301.52019326729472951, ERR);

        // (3,1,0,0,0,...)
        assert_approx_eq!(orbitBA(gen.get_l3(), gen.get_l1(), &cube5d, f_1234 , 5, &mut f_evals), 556.36420068745234470, ERR);

        // (3,2,0,0,0,...)
        assert_approx_eq!(orbitBA(gen.get_l3(), gen.get_l2(), &cube5d, f_1234 , 5, &mut f_evals), 275.91878544684016599, ERR);

        // (4,1,0,0,0,...)
        assert_approx_eq!(orbitBA(gen.get_l4(), gen.get_l1(), &cube5d, f_1234 , 5, &mut f_evals), 521.62779899524017755, ERR);


    }

    #[test]
    fn orbit_BAA() {
        let mut f_evals:usize = 0;
        // Initialize generators
        let gen = Generators::new(100);
        let cube5d = HyperCube::new(&vec![1.0; 5], &vec![1.0; 5]);
        const ERR:f64 = 1e-9f64;

        // (2,1,1,0,0,...)
        assert_approx_eq!(orbitBAA(gen.get_l2(), gen.get_l1(), &cube5d, f_1234 , 5, &mut f_evals), 2169.0227707260377438, ERR);

        // (3,1,1,0,0,...)
        assert_approx_eq!(orbitBAA(gen.get_l3(), gen.get_l1(), &cube5d, f_1234 , 5, &mut f_evals), 3628.2529347296686588, ERR);
    }


    #[test]
    fn orbit_BBA() {            // Small hack!!! A <-> B
        let mut f_evals:usize = 0;
        // Initialize generators
        let gen = Generators::new(100);
        let cube5d = HyperCube::new(&vec![1.0; 5], &vec![1.0; 5]);
        const ERR:f64 = 1e-9f64;

        // (1,2,2,0,0,...)
        assert_approx_eq!(orbitBAA(gen.get_l1(), gen.get_l2(), &cube5d, f_1234 , 5, &mut f_evals), 1119.9590592233820280, ERR);

        assert_approx_eq!(orbitBAA(gen.get_l1(), gen.get_l2(), &cube5d, f_2222 , 5, &mut f_evals), 528.44703227114887674, ERR);

    }

    #[test]
    fn orbit_BAAA() {
        let mut f_evals:usize = 0;
        // Initialize generators
        let gen = Generators::new(100);
        let cube5d = HyperCube::new(&vec![1.0; 5], &vec![1.0; 5]);
        const ERR:f64 = 1e-9f64;

        // (2,1,1,0,0,...)
        assert_approx_eq!(orbitBAAA(gen.get_l2(), gen.get_l1(), &cube5d, f_1234 , 5, &mut f_evals), 5991.4976627978102740, ERR);

    }

}
