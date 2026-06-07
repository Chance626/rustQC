use faer::Mat;

pub fn print_2D_mat(twodmat: &Mat::<f64>) {
    let print_twodmat: Vec<Vec<f64>> = (0..twodmat.nrows())
        .map(|r| (0..twodmat.ncols()).map(|c| twodmat[(r, c)]).collect()).collect();

    for row in &print_twodmat {
        for i in row {
            print!("{:<0.4} ", i);
        }
        println!();
    }

}