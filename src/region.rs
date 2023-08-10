#![allow(non_snake_case)]
// use crate::integrator::Point;
use crate::integrator::IntegrandFunction;
use crate::integrator::HyperCube;
use crate::integrator::RegionResult;
use crate::rules::Generators;
use crate::rules::Rule11weights;
use crate::rules::Rule13weights;
// Orbits
use crate::orbits as orb;

pub fn integrate_region(hc:&HyperCube, f:IntegrandFunction, n:usize, r11:&Rule11weights, r13:&Rule13weights, g:&Generators, _id:&'static str) -> RegionResult {
    // perform function evaluation across different orbits
    let mut f_evals:usize = 0;

    // vectors to store functions with steps along axes
    let mut f_lam1:Vec<(f64,f64)> = vec![(0f64,0f64); n];
    let mut f_lam2:Vec<(f64,f64)> = vec![(0f64,0f64); n];
    let mut f_lam3:Vec<(f64,f64)> = vec![(0f64,0f64); n];
    let mut f_lam4:Vec<(f64,f64)> = vec![(0f64,0f64); n];
    let mut f_lam5:Vec<(f64,f64)> = vec![(0f64,0f64); n];

    let mut res11:f64 = 0f64;
    let mut res13:f64 = 0f64;

    // (0,0,0,0,0,...)
    let f0 = f(&hc.xc);
    f_evals +=1;

    let orb_00000 = f0*hc.hwvol;
    res11 += r11.get_w00000()*orb_00000;
    res13 += r13.get_w00000()*orb_00000;

    // (1,0,0,0,0,...)
    // let orb_10000 = orb::orbitAdiff(g.get_l1(), &hc, f , n, &mut lam1p, &mut lam1m);
    let orb_10000 = orb::orbitAdiff(g.get_l1(), hc, f , n, &mut f_lam1, &mut f_evals);
    res11 += r11.get_w10000()*orb_10000;
    res13 += r13.get_w10000()*orb_10000;

    // (1,1,0,0,0,...)
    let orb_11000 = orb::orbitAA(g.get_l1(), hc, f , n, &mut f_evals);
    res11 += r11.get_w11000()*orb_11000;
    res13 += r13.get_w11000()*orb_11000;

    // (1,1,1,0,0,...)
    let orb_11100 = orb::orbitAAA(g.get_l1(), hc, f , n, &mut f_evals);
    res11 += r11.get_w11100()*orb_11100;
    res13 += r13.get_w11100()*orb_11100;

    // (1,1,1,1,0,...)
    let orb_11110 = orb::orbitAAAA(g.get_l1(), hc, f , n, &mut f_evals);
    res11 += r11.get_w11110()*orb_11110;
    res13 += r13.get_w11110()*orb_11110;

    // (1,1,1,1,1,...)
    let orb_11111 = orb::orbitAAAAA(g.get_l1(), hc, f , n, &mut f_evals);
    // res11 += r11.get_w11111()*orb_11111;
    res13 += r13.get_w11111()*orb_11111;

    // (2,0,0,0,0,...)
    // let orb_20000 = orb::orbitAdiff(g.get_l2(), hc, f , n, &mut lam2p, &mut lam2m);
    let orb_20000 = orb::orbitAdiff(g.get_l2(), hc, f , n, &mut f_lam2, &mut f_evals);
    res11 += r11.get_w20000()*orb_20000;
    res13 += r13.get_w20000()*orb_20000;

    // (2,1,0,0,0,...)
    let orb_21000 = orb::orbitBA(g.get_l2(), g.get_l1(), hc, f , n, &mut f_evals);
    res11 += r11.get_w21000()*orb_21000;
    res13 += r13.get_w21000()*orb_21000;

    // (2,1,1,0,0,...)
    let orb_21100 = orb::orbitBAA(g.get_l2(), g.get_l1(), hc, f , n, &mut f_evals);
    res11 += r11.get_w21100()*orb_21100;
    res13 += r13.get_w21100()*orb_21100;

    // (2,1,1,1,0,...)
    let orb_21110 = orb::orbitBAAA(g.get_l2(), g.get_l1(), hc, f , n, &mut f_evals);
    // res11 += r11.get_w21110()*orb_21110;
    res13 += r13.get_w21110()*orb_21110;

    // (2,2,0,0,0,...)
    // res += r.get_w22000()*
    let orb_22000 = orb::orbitAA(g.get_l2(), hc, f , n, &mut f_evals);
    res11 += r11.get_w22000()*orb_22000;
    res13 += r13.get_w22000()*orb_22000;

    // (2,2,1,0,0,...) TODO
    let orb_22100 = orb::orbitBAA(g.get_l1(), g.get_l2(), hc, f , n, &mut f_evals);
    // res11 += r11.get_w22100()*orb_22100;
    res13 += r13.get_w22100()*orb_22100;

    // (3,0,0,0,0,...)
    // let orb_30000 = orb::orbitAdiff(g.get_l3(), hc, f , n, &mut lam3p, &mut lam3m);
    let orb_30000 = orb::orbitAdiff(g.get_l3(), hc, f , n, &mut f_lam3, &mut f_evals);
    res11 += r11.get_w30000()*orb_30000;
    res13 += r13.get_w30000()*orb_30000;

    // (3,1,0,0,0,...)
    let orb_31000 = orb::orbitBA(g.get_l3(), g.get_l1(), hc, f , n, &mut f_evals);
    res11 += r11.get_w31000()*orb_31000;
    res13 += r13.get_w31000()*orb_31000;

    // (3,1,1,0,0,...)
    let orb_31100 = orb::orbitBAA(g.get_l3(), g.get_l1(), hc, f , n, &mut f_evals);
    // res11 += r11.get_w31100()*orb_31100;
    res13 += r13.get_w31100()*orb_31100;

    // (3,2,0,0,0,...)
    let orb_32000 = orb::orbitBA(g.get_l3(), g.get_l2(), hc, f , n, &mut f_evals);
    // res11 += r11.get_w32000()*orb_32000;
    res13 += r13.get_w32000()*orb_32000;

    // (4,0,0,0,0,...)
    // let orb_40000 = orb::orbitAdiff(g.get_l4(), hc, f , n, &mut lam4p, &mut lam4m);
    let orb_40000 = orb::orbitAdiff(g.get_l4(), hc, f , n, &mut f_lam4, &mut f_evals);
    res11 += r11.get_w40000()*orb_40000;
    res13 += r13.get_w40000()*orb_40000;

    // (4,1,0,0,0,...)
    let orb_41000 = orb::orbitBA(g.get_l4(), g.get_l1(), hc, f , n, &mut f_evals);
    // res11 += r11.get_41000()*orb_41000;
    res13 += r13.get_w41000()*orb_41000;

    // (5,0,0,0,0,...)
    // let orb_50000 = orb::orbitAdiff(g.get_l5(), hc, f , n, &mut lam5p, &mut lam5m);
    let orb_50000 = orb::orbitAdiff(g.get_l5(), hc, f , n, &mut f_lam5, &mut f_evals);
    // res11 += r11.get_w50000()*orb_50000;
    res13 += r13.get_w50000()*orb_50000;

    // (delta,delta,...,delta)
    // res += r.get_wDelta()*
    let orb_delta = orb::orbitDelta(g.get_de(), hc, f , n, &mut f_evals);
    res11 += r11.get_w_delta()*orb_delta;
    res13 += r13.get_w_delta()*orb_delta;

    // Classical fourth difference
    let ratio = g.get_l1s()/g.get_l2s();
    let mut maxD:f64   = 0.0;
    let mut idxD:usize = 0;
    // println!("l1/l2 r = {}",ratio);
    for (i, f_a1_i) in f_lam1.iter().enumerate() {
        let di_f = (f_a1_i.0 + f_a1_i.1 - 2.0*f0 - ratio*(f_lam2[i].0 + f_lam2[i].1 - 2.0*f0)).abs();
        if di_f > maxD {maxD = di_f; idxD = i}
        // println!("Element at position {}: {}, info:: f0={} f1p={} f1m={} ", i, di_f, f0, f_a1_i.0,f_a1_i.1);
    }
    // println!("max(D) = {}, i = {}", maxD,idxD);

    // // Classical fourth difference
    // let ratio = g.get_l1s()/g.get_l3s();
    // println!("l1/l3 r = {}",ratio);
    // for (i, f_a1_i) in f_lam1.iter().enumerate() {
    //     let di_f = (f_a1_i.0 + f_a1_i.1 - 2.0*f0 - ratio*(f_lam3[i].0 + f_lam3[i].1 - 2.0*f0)).abs();
    //     println!("Element at position {}: {:?}", i, di_f);
    // }


    // println!("ERR = {}", (res11/res13-1f64).abs());
    // (res11, res13, idxD)
    RegionResult {
        region: hc.clone(),
        res11,
        res13,
        res: res13,
        err: (res11-res13).abs(),
        divdir: idxD,
        f_evals,
    }
}




pub fn split_region(r0:&HyperCube, k:usize) -> (HyperCube, HyperCube) {
    let mut r_l = r0.clone();
    let mut r_r = r0.clone();

    r_l.xc[k] -= r_l.hw[k]/2f64;
    r_l.hw[k] *= 0.5f64;
    r_l.hwvol *= 0.5f64;

    r_r.xc[k] += r_r.hw[k]/2f64;
    r_r.hw[k] *= 0.5f64;
    r_r.hwvol *= 0.5f64;

    // TODO add to lng list for calulation
    (r_l,r_r)
}

#[cfg(test)]
mod tests {
    use crate::integrator::Point;
    fn int_rule13(hc:&HyperCube, f:IntegrandFunction, n:usize, r:&Rule13weights, g:&Generators) -> f64 {
        let mut f_evals:usize = 0;
        // (0,0,0,0,0,...)
        let mut res:f64 = r.get_w00000()*f(&hc.xc)*hc.hwvol;
        // (1,0,0,0,0,...)
        res += r.get_w10000()*orb::orbitA(g.get_l1(), hc, f , n , &mut f_evals);
        // (1,1,0,0,0,...)
        res += r.get_w11000()*orb::orbitAA(g.get_l1(), hc, f , n , &mut f_evals);
        // (1,1,1,0,0,...)
        res += r.get_w11100()*orb::orbitAAA(g.get_l1(), hc, f , n , &mut f_evals);
        // (1,1,1,1,0,...)
        res += r.get_w11110()*orb::orbitAAAA(g.get_l1(), hc, f , n , &mut f_evals);
        // (1,1,1,1,1,...)
        res += r.get_w11111()*orb::orbitAAAAA(g.get_l1(), hc, f , n , &mut f_evals);
        // (2,0,0,0,0,...)
        res += r.get_w20000()*orb::orbitA(g.get_l2(), hc, f , n , &mut f_evals);
        // (2,1,0,0,0,...)
        res += r.get_w21000()*orb::orbitBA(g.get_l2(), g.get_l1(), hc, f , n , &mut f_evals);
        // (2,1,1,0,0,...)
        res += r.get_w21100()*orb::orbitBAA(g.get_l2(), g.get_l1(), hc, f , n , &mut f_evals);
        // (2,1,1,1,0,...)
        res += r.get_w21110()*orb::orbitBAAA(g.get_l2(), g.get_l1(), hc, f , n , &mut f_evals);
        // (2,2,0,0,0,...)
        res += r.get_w22000()*orb::orbitAA(g.get_l2(), hc, f , n , &mut f_evals);
        // (2,2,1,0,0,...) TODO
        res += r.get_w22100()*orb::orbitBAA(g.get_l1(), g.get_l2(), hc, f , n , &mut f_evals);
        // (3,0,0,0,0,...)
        res += r.get_w30000()*orb::orbitA(g.get_l3(), hc, f , n , &mut f_evals);
        // (3,1,0,0,0,...)
        res += r.get_w31000()*orb::orbitBA(g.get_l3(), g.get_l1(), hc, f , n , &mut f_evals);
        // (3,1,1,0,0,...)
        res += r.get_w31100()*orb::orbitBAA(g.get_l3(), g.get_l1(), hc, f , n , &mut f_evals);
        // (3,2,0,0,0,...)
        res += r.get_w32000()*orb::orbitBA(g.get_l3(), g.get_l2(), hc, f , n , &mut f_evals);
        // (4,0,0,0,0,...)
        res += r.get_w40000()*orb::orbitA(g.get_l4(), hc, f , n , &mut f_evals);
        // (4,1,0,0,0,...)
        res += r.get_w41000()*orb::orbitBA(g.get_l4(), g.get_l1(), hc, f , n , &mut f_evals);
        // (5,0,0,0,0,...)
        res += r.get_w50000()*orb::orbitA(g.get_l5(), hc, f , n , &mut f_evals);
        // (delta,delta,...,delta)
        res += r.get_w_delta()*orb::orbitDelta(g.get_de(), hc, f , n , &mut f_evals);

        res
    }


    fn int_rule11(hc:&HyperCube, f:IntegrandFunction, n:usize, r:&Rule11weights, g:&Generators) -> f64 {
        let mut f_evals:usize = 0;
        // (0,0,0,0,0,...)
        let mut res:f64 = r.get_w00000()*f(&hc.xc)*hc.hwvol;
        // (1,0,0,0,0,...)
        res += r.get_w10000()*orb::orbitA(g.get_l1(), hc, f , n , &mut f_evals );
        // (1,1,0,0,0,...)
        res += r.get_w11000()*orb::orbitAA(g.get_l1(), hc, f , n , &mut f_evals );
        // (1,1,1,0,0,...)
        res += r.get_w11100()*orb::orbitAAA(g.get_l1(), hc, f , n , &mut f_evals );
        // (1,1,1,1,0,...)
        res += r.get_w11110()*orb::orbitAAAA(g.get_l1(), hc, f , n , &mut f_evals );
        // (2,0,0,0,0,...)
        res += r.get_w20000()*orb::orbitA(g.get_l2(), hc, f , n , &mut f_evals );
        // (2,1,0,0,0,...)
        res += r.get_w21000()*orb::orbitBA(g.get_l2(), g.get_l1(), hc, f , n , &mut f_evals );
        // (2,1,1,0,0,...)
        res += r.get_w21100()*orb::orbitBAA(g.get_l2(), g.get_l1(), hc, f , n , &mut f_evals );
        // (2,2,0,0,0,...)
        res += r.get_w22000()*orb::orbitAA(g.get_l2(), hc, f , n , &mut f_evals );
        // (3,0,0,0,0,...)
        res += r.get_w30000()*orb::orbitA(g.get_l3(), hc, f , n , &mut f_evals );
        // (3,1,0,0,0,...)
        res += r.get_w31000()*orb::orbitBA(g.get_l3(), g.get_l1(), hc, f , n , &mut f_evals );
        // (4,0,0,0,0,...)
        res += r.get_w40000()*orb::orbitA(g.get_l4(), hc, f , n , &mut f_evals );
        // (delta,delta,...,delta)
        res += r.get_w_delta()*orb::orbitDelta(g.get_de(), hc, f , n , &mut f_evals );

        res
    }




    use assert_approx_eq::assert_approx_eq;
    use super::*;

    // Test functions
    fn f_0(_:&Point) -> f64 { 1.0 }
    fn f_11(x:&Point) -> f64 { x[0]*x[1] }
    fn f_12(x:&Point) -> f64 { x[0]*x[1].powi(2) }
    fn f_22(x:&Point) -> f64 { x[0].powi(2)*x[1].powi(2) }
    fn f_222(x:&Point) -> f64 { x[0].powi(2)*x[1].powi(2)*x[2].powi(2) }
    fn f_2222(x:&Point) -> f64 { x[0].powi(2)*x[1].powi(2)*x[2].powi(2)*x[3].powi(2) }
    fn f_1234(x:&Point) -> f64 { x[0].powi(1)*x[1].powi(2)*x[2].powi(3)*x[3].powi(4) }
    fn f_23232(x:&Point) -> f64 { x[0].powi(2)*x[1].powi(3)*x[2].powi(2)*x[3].powi(3)*x[4].powi(2) }

    const ERR:f64 = 1e-9f64;

    #[test]
    fn rule11_dim5() {

        const DIM:usize = 5;
        let r11 = Rule11weights::new(100, DIM);
        let gen = Generators::new(100);
        // [0,2]
        let cube5d02 = HyperCube::new(&vec![1.0; DIM], &vec![1.0; DIM]);
        // [-1,1]
        let cube5dm1p1 = HyperCube::new(&vec![0.0; DIM], &vec![1.0; DIM]);
        // [-2,2]
        let cube5dm2p2 = HyperCube::new(&vec![0.0; DIM], &vec![2.0; DIM]);
        // [0,1]
        let cube5d01 = HyperCube::new(&vec![0.5; DIM], &vec![0.5; DIM]);

        // Standard cube
        assert_approx_eq!(int_rule11(&cube5dm1p1, f_22, DIM, &r11, &gen), 32f64/9f64, ERR);
        // Shifted cube
        assert_approx_eq!(int_rule11(&cube5d02, f_22, DIM, &r11, &gen), 512f64/9f64, ERR);

        // Rescaled cube
        assert_approx_eq!(int_rule11(&cube5dm2p2, f_0, DIM, &r11, &gen), 1024f64, ERR);
        assert_approx_eq!(int_rule11(&cube5dm2p2, f_22, DIM, &r11, &gen), 16384f64/9f64, ERR);

        // Real life [0,1] cube
        assert_approx_eq!(int_rule11(&cube5d01, f_22, DIM, &r11, &gen), 1f64/9f64, ERR);
        assert_approx_eq!(int_rule11(&cube5d01, f_0, DIM, &r11, &gen), 1f64, ERR);
        assert_approx_eq!(int_rule11(&cube5d01, f_1234, DIM, &r11, &gen), 1f64/120f64, ERR);
    }

    #[test]
    fn rule13_dim5() {

        const DIM:usize = 5;
        let r13 = Rule13weights::new(100, DIM);
        let gen = Generators::new(100);
        // [0,2]
        let cube5d02 = HyperCube::new(&vec![1.0; DIM], &vec![1.0; DIM]);
        // [-1,1]
        let cube5dm1p1 = HyperCube::new(&vec![0.0; DIM], &vec![1.0; DIM]);
        // [-2,2]
        let cube5dm2p2 = HyperCube::new(&vec![0.0; DIM], &vec![2.0; DIM]);
        // [0,1]
        let cube5d01 = HyperCube::new(&vec![0.5; DIM], &vec![0.5; DIM]);

        // Standard cube
        assert_approx_eq!(int_rule13(&cube5dm1p1, f_22, DIM, &r13, &gen), 32f64/9f64, ERR);
        // Shifted cube
        assert_approx_eq!(int_rule13(&cube5d02, f_22, DIM, &r13, &gen), 512f64/9f64, ERR);

        // Rescaled cube
        assert_approx_eq!(int_rule13(&cube5dm2p2, f_0, DIM, &r13, &gen), 1024f64, ERR);
        assert_approx_eq!(int_rule13(&cube5dm2p2, f_22, DIM, &r13, &gen), 16384f64/9f64, ERR);

        // Real life [0,1] cube
        assert_approx_eq!(int_rule13(&cube5d01, f_22, DIM, &r13, &gen), 1f64/9f64, ERR);
        assert_approx_eq!(int_rule13(&cube5d01, f_0, DIM, &r13, &gen), 1f64, ERR);
        assert_approx_eq!(int_rule13(&cube5d01, f_1234, DIM, &r13, &gen), 1f64/120f64, ERR);
    }


    #[test]
    fn integrate_region_dim5() {

        const DIM:usize = 5;
        let r11 = Rule11weights::new(100, DIM);
        let r13 = Rule13weights::new(100, DIM);
        let gen = Generators::new(100);
        // [0,2]
        let _cube5d02 = HyperCube::new(&vec![1.0; DIM], &vec![1.0; DIM]);
        // [-1,1]
        let _cube5dm1p1 = HyperCube::new(&vec![0.0; DIM], &vec![1.0; DIM]);
        // [-2,2]
        let _cube5dm2p2 = HyperCube::new(&vec![0.0; DIM], &vec![2.0; DIM]);
        // [0,1]
        let cube5d01 = HyperCube::new(&vec![0.5; DIM], &vec![0.5; DIM]);

        // Standard cube
        // assert_approx_eq!(integrate_region(&cube5dm1p1, f_22, DIM, &r11, &r13, &gen), 32f64/9f64, ERR);
        // // Shifted cube
        // assert_approx_eq!(int_rule11(&cube5d02, f_22, DIM, &r11, &gen), 512f64/9f64, ERR);

        // // Rescaled cube
        // assert_approx_eq!(int_rule11(&cube5dm2p2, f_0, DIM, &r11, &gen), 1024f64, ERR);
        // assert_approx_eq!(int_rule11(&cube5dm2p2, f_22, DIM, &r11, &gen), 16384f64/9f64, ERR);

        // // Real life [0,1] cube
        // assert_approx_eq!(int_rule11(&cube5d01, f_22, DIM, &r11, &gen), 1f64/9f64, ERR);
        // assert_approx_eq!(int_rule11(&cube5d01, f_0, DIM, &r11, &gen), 1f64, ERR);
        // assert_approx_eq!(integrate_region(&cube5d01, f_1234, DIM, &r11, &r13, &gen), 1f64/120f64, ERR);


        assert_approx_eq!(integrate_region(&cube5d01, f_1234, DIM, &r11, &r13, &gen,"").res11, 1f64/120f64, ERR);

        // we need a new cube since integrate_region destroys existing one
        // let cube5d01 = HyperCube::new(&vec![0.5; DIM], &vec![0.5; DIM]);
        assert_approx_eq!(integrate_region(&cube5d01, f_1234, DIM, &r11, &r13, &gen,"").res13, 1f64/120f64, ERR);


        let _cube5d01 = HyperCube::new(&vec![0.5; DIM], &vec![0.5; DIM]);

    }


}
