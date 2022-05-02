

fn f_30ord(x:&incubator::Point) -> f64 { (x[0].powi(1)*x[1].powi(2)*x[2].powi(3)*x[3].powi(4)*x[4].powi(5)).powi(2) }


fn main() {


    let mut ii = incubator::Integrator::new(f_30ord, &vec![0.0; 5], &vec![1.0; 5]);

    let num_cores = 6;

    ii.calc(num_cores,
            2000,                 // max regions
            100000,                // max f calls
            0.0001                // rel err
    );

    ii.dump();
}
