

fn f_15ord(x:&incubator::Point) -> f64 { (x[0]+1.0)/(x[0]*x[1]-1.0)*x[1].powi(2)*x[2].powi(3)*x[3].powi(4)*x[4].powi(5) }

fn main() {


    let mut ii = incubator::Integrator::new(f_15ord, &vec![0.0; 5], &vec![1.0; 5]);

    let num_cores = 6;

    ii.calc(num_cores,
            2000,                 // max regions
            100000,                // max f calls
            0.0001                // rel err
    );

    ii.dump();
}
