use crate::rules;
use crate::region::integrate_region;
use crate::region::split_region;


use std::collections::BinaryHeap;
use std::cmp::Ordering;
use std::cmp;

use std::collections::HashMap;

use std::sync::atomic::AtomicUsize;
use std::sync::atomic::Ordering as atmord; // to fix name overlapping

use std::time::Instant;

use serde::{Serialize, Deserialize};
// use serde_with::serde_as;

// Reasons to break integration loop
pub enum IntegrationStopReason {
    RelErr,                     // Relative error goal achieved
    AbsErr,                     // Absolute
    NumFunCalls,                // Too many function calls
    NumRegions,                 // Too many regions

    InternalError,              // Error
}


const NREGIONSALLOC:usize = 1000;

// n-dimensional vector
pub type Point = Vec<f64>;
// type Point = [Float];
// function signature to be integrated
pub type IntegrandFunction = fn(&Point) -> f64;

#[derive(Clone, Serialize, Deserialize, Debug)]
pub struct HyperCube {
    pub xc: Point,                  // Center coordinate
    pub hw: Point,                  // Vector of half-widths
    pub hwvol: f64,                 // Volume from half-widths
    pub name: String,
}


impl Eq for HyperCube {}
impl PartialEq for HyperCube {
    fn eq(&self, _other: &Self) -> bool {false}
}

impl HyperCube {
    pub fn new(xc:&Point, hw:&Point) -> HyperCube {
        let mut vol:f64 = 1.0;
        for w in hw.iter() { vol *= w};
        HyperCube {
            xc:xc.to_vec(),
            hw:hw.to_vec(),
            hwvol:vol,
            name:"".to_string(),
        }
    }

    pub fn from_bounds(lower_bounds:&Point, upper_bounds:&Point) -> HyperCube {
        assert!(lower_bounds.len() == upper_bounds.len(),
                "dim(a) = {}, dim(b) = {}", lower_bounds.len(), upper_bounds.len());
        let mut center  = Point::new();
        let mut widths = Point::new();
        let mut vol:f64 = 1f64;
        for (&a_i, &b_i) in lower_bounds.iter().zip(upper_bounds.iter()) {
            center.push((a_i+b_i)/2f64);
            widths.push((b_i-a_i)/2f64);
            vol *= (b_i-a_i)/2f64;
            // println!("a = {}, b = {}", a_i, b_i);
        }
        HyperCube {
            xc:center.to_vec(),
            hw:widths.to_vec(),
            hwvol:vol,
            name:"".to_string(),
        }
    }

    pub fn set_name(&mut self, name:&String) -> &HyperCube {
        self.name = name.to_string();
        self
    }

}


#[derive(Serialize, Deserialize, Debug)]
pub struct RegionResult {
    pub region: HyperCube,
    pub res11: f64,                 // result from Rule11
    pub res13: f64,                 // result from Rule13
    pub res: f64,                   // final result
    pub err: f64,                   // error for further subdivision
    pub divdir: usize,              // direction of the subdivision
    pub f_evals: usize,             // number of function evaluations in region
}

impl Eq for RegionResult {}

impl Ord for RegionResult {
    fn cmp(&self, other: &Self) -> Ordering {
        self.err.partial_cmp(&other.err).unwrap() // FIXME check for NaN compare
    }
}

impl PartialOrd for RegionResult {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        self.err.partial_cmp(&other.err)
    }
}

impl PartialEq for RegionResult {
    fn eq(&self, other: &Self) -> bool {
        self.err == other.err
    }
}


pub struct Integrator {
    ndim: usize,
    f: IntegrandFunction,
    start_region: HyperCube,
    // Heap of calculated regions, sorted by error in descending order
    // pop brings region with the highest error
    regions_heap: BinaryHeap<RegionResult>,
    r11: rules::Rule11weights,
    r13: rules::Rule13weights,
    gen: rules::Generators,
    // Format: [reg_id, (reg_res, left_id, right_id)]
    saved_regions: HashMap<String,(RegionResult,String,String)>,
    // final parameters
    f_calls: usize,       // number of function calls
    n_regions: usize,     // number of regions after subdivision
    res: f64,             // final result
    abserr: f64,          // absolute error
    relerr: f64,          // relative error
}

impl Integrator {
    pub fn new(f:IntegrandFunction, lower_bounds:&Point, upper_bounds:&Point)  -> Integrator{
        assert!(lower_bounds.len() == upper_bounds.len(),
                "dim(a) = {}, dim(b) = {}", lower_bounds.len(), upper_bounds.len());
        let dim = lower_bounds.len();
        Integrator {
            f,
            ndim: dim,
            start_region: HyperCube::from_bounds(lower_bounds, upper_bounds),
            regions_heap: BinaryHeap::with_capacity(NREGIONSALLOC),
            r11: rules::Rule11weights::new(100, dim),
            r13: rules::Rule13weights::new(100, dim),
            gen: rules::Generators::new(100),
            saved_regions: HashMap::new(),
            // final parameters
            f_calls: 0,       // number
            n_regions: 0,
            abserr: 0.0,
            relerr: 0.0,
            res: 0.0,
        }
    }


    pub fn calc(&mut self, ncores:usize, max_regions:usize, max_f_calls:usize, rel_err_goal:f64) {
        println!("Integration started with ncores = {ncores}");
        let start_time = Instant::now();
        // Store number of function calls
        let f_evals = AtomicUsize::new(0);

        rayon::ThreadPoolBuilder::new().num_threads(ncores).build_global().unwrap();

        // We make a copy of the starting region to keep track of the integration region bounds
        // let start_cube = self.start_region.clone();
        // make first integration and start using a heap
        let sr:RegionResult = integrate_region(self.start_region.set_name(&"S_".to_string()), self.f, self.ndim, &self.r11, &self.r13, &self.gen, "555");
        f_evals.fetch_add(sr.f_evals, atmord::SeqCst);
        self.relerr = (sr.err/sr.res).abs();
        self.abserr = sr.err;
        self.res = sr.res;
        self.regions_heap.push(sr);

        // Main loop for regions subdivision
        // Termination conditions:
        //   1. max number of regions exceeded
        //   2. number of function evaluations limit reached
        //   3. relative error achieved
        //   4. absolute error achieved
        while (self.regions_heap.len() < max_regions) && ( self.f_calls < max_f_calls) && (self.relerr > rel_err_goal) {
            // TODO better preallocation
            let mut new_regions: Vec<HyperCube> = Vec::with_capacity(ncores);

            for _ in 0..cmp::min(ncores/2, self.regions_heap.len())
            {
                // println!("Popping region!");
                let r = self.regions_heap.pop().unwrap();

                let dir = r.divdir;
                let killed_region_name = r.region.name.clone();
                // Split calculated region cube in two left(L) and right(R)
                let (mut cube_l, mut cube_r) = split_region(&r.region, dir);
                let region_l_name = format!("{killed_region_name}{dir}L");
                let region_r_name = format!("{killed_region_name}{dir}R");

                cube_l.set_name(&region_l_name);
                cube_r.set_name(&region_r_name);

                new_regions.push(cube_l);
                new_regions.push(cube_r);

                // Save region result: [reg_id, (reg_res, left_id, right_id)]
                self.saved_regions.insert(killed_region_name, (r, region_l_name, region_r_name));
            }

            use rayon::prelude::*;
            let mut result: BinaryHeap<RegionResult> = BinaryHeap::with_capacity(ncores);

            result = new_regions.par_iter().map(|cube| {
                // calculation result to save f_calls
                let rr:RegionResult = integrate_region(cube, self.f, self.ndim, &self.r11, &self.r13, &self.gen,"666");
                f_evals.fetch_add(rr.f_evals, atmord::SeqCst);
                rr
                }).collect();
            self.regions_heap.append(&mut result);
            // println!("new regions length = {}", self.regions_heap.len());


            // Update total number of function calls
            self.f_calls = f_evals.load(atmord::SeqCst);
            // Update number of regions
            self.n_regions = self.regions_heap.len();      // TODO
            // Update error
            let mut total:f64 = 0.0;
            let mut err2 :f64 = 0.0;
            for x in self.regions_heap.iter() {
                total += x.res;
                err2  += (x.res11-x.res13)*(x.res11-x.res13);
            }
            self.abserr = err2.sqrt();
            self.relerr = (err2.sqrt()/total).abs();
            self.res    = total;

            println!(" *  N(regions)={regions:>6}  || log10(fcalls)={calls:>10.4} || sigma={sigma:10.4}  || I=({res:+12.4e} +/- {err:12.4e})",regions=self.n_regions, calls=(self.f_calls as f64).log10(), sigma = self.relerr,res=self.res, err = self.abserr);
            println!();
            println!("    Elapsed time: {:.2?}", start_time.elapsed());
            println!();
        }

        let mut stop_reason:IntegrationStopReason = IntegrationStopReason::InternalError;

        if self.regions_heap.len() >= max_regions {
            stop_reason = IntegrationStopReason::NumRegions;
        } else if self.f_calls >= max_f_calls {
            stop_reason = IntegrationStopReason::NumFunCalls;
        } else if self.relerr <= rel_err_goal {
            stop_reason = IntegrationStopReason::RelErr;
        }

        // Final statistics
        println!();
        println!();
        println!();
        println!("===========================================================================");
        match stop_reason {
            IntegrationStopReason::RelErr        => println!("stop reason       = RelErr"),
            IntegrationStopReason::AbsErr        => println!("stop reason       = AbsErr"),
            IntegrationStopReason::NumFunCalls   => println!("stop reason       = NumFunCalls"),
            IntegrationStopReason::NumRegions    => println!("stop reason       = NumRegions"),
            IntegrationStopReason::InternalError => println!("stop reason       = ERROR"),
        }
        println!("function calls    = {}", self.f_calls);
        println!("number of regions = {}", self.n_regions);
        println!("relative error    = {}", self.relerr);
        println!("total result      = {} +/- {}", self.res, self.abserr);
        println!("===========================================================================");
    }



    pub fn dump(&self) {

        // TODO
    }

}






// -------------------------------------------------------------------------------------------------------------
//
//                                            TESTS PART
//
// -------------------------------------------------------------------------------------------------------------
#[cfg(test)]
mod tests {

    use assert_approx_eq::assert_approx_eq;
    use super::*;

    const ERR:f64 = 1e-9f64;

    fn f_15ord(x:&Point) -> f64 {


        (x[0]+1.0)/(x[0]*x[1]-1.0)*x[1].powi(2)*x[2].powi(3)*x[3].powi(4)*x[4].powi(5) }

    fn f_30ord(x:&Point) -> f64 { (x[0].powi(1)*x[1].powi(2)*x[2].powi(3)*x[3].powi(4)*x[4].powi(5)).powi(2) }

    #[test]
    fn hcube_bounds() {

        let hc = HyperCube::from_bounds(&vec![0.0; 3], &vec![13.0; 3]);

        assert_approx_eq!(hc.hwvol, 2197f64/8f64, ERR);
    }

    #[test]
    fn integrator_init() {

        let mut ii = Integrator::new(f_15ord, &vec![0.0; 5], &vec![1.0; 5]);

        ii.calc(3,200,1000,0.001);

        ii.dump();
    }

}
