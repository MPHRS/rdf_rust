#![deny(clippy::all)]
use plotters::prelude::*;
use std::error::Error;
use std::fs::File;
use std::io::Write;
use std::io::{BufRead, BufReader};
struct MyBox {
    x: f32,
    y: f32,
    z: f32,
}

impl MyBox {
    fn new(x: f32, y: f32, z: f32) -> MyBox {
        MyBox { x, y, z }
    }

    fn volume(&self) -> f32 {
        self.x * self.y * self.z
    }

    fn periodic_correct(&self, xb: f32, yb: f32, zb: f32) -> (f32, f32, f32) {
        let xb = MyBox::periodic(xb, self.x);
        let yb = MyBox::periodic(yb, self.y);
        let zb = MyBox::periodic(zb, self.z);
        (xb, yb, zb)
    }

    fn check_in_box(&self, xb: f32, yb: f32, zb: f32) -> bool {
        (xb.abs() < 0.5 * self.x) && (yb.abs() < 0.5 * self.y) && (zb.abs() < 0.5 * self.z)
    }

    fn periodic(coord: f32, box_size: f32) -> f32 {
        if coord.abs() > 1.5 * box_size {
            panic!("Out of box");
        }
        if coord.abs() > 0.5 * box_size {
            coord - coord.signum() * box_size
        } else {
            coord
        }
    }
}
struct Atom {
    x: f32,
    y: f32,
    z: f32,
    btype: i32,
}

impl Atom {
    fn new(x: f32, y: f32, z: f32, btype: i32) -> Atom {
        Atom { x, y, z, btype }
    }

    fn get_distance(&self, atom: &Atom, box_obj: &MyBox) -> f32 {
        let dx = atom.x - self.x;
        let dy = atom.y - self.y;
        let dz = atom.z - self.z;
        let (dx, dy, dz) = box_obj.periodic_correct(dx, dy, dz);
        (dx.powi(2) + dy.powi(2) + dz.powi(2)).sqrt()
    }
}

struct ReadTrack {
    path: String,
    time_step: i32,
    file_track: BufReader<File>,
    num_atoms: i32,
    my_box: MyBox,
    x: Vec<f32>,
    y: Vec<f32>,
    z: Vec<f32>,
    btype: Vec<i32>,
    iterator: Vec<usize>,
}

impl ReadTrack {
    fn new(path: &str) -> Result<ReadTrack, std::io::Error> {
        let full_path = format!("{}/{}", path, "TRACK");
        let file = File::open(&full_path)?;
        let mut title = String::new();
        let mut file_track = BufReader::new(file);
        file_track.read_line(&mut title)?;
        let title_parts: Vec<&str> = title.split_whitespace().collect();
        let num_atoms = title_parts[1]
            .parse::<i32>()
            .map_err(|_| std::io::ErrorKind::InvalidData)?;
        let my_box = MyBox::new(
            title_parts[3]
                .parse::<f32>()
                .map_err(|_| std::io::ErrorKind::InvalidData)?,
            title_parts[4]
                .parse::<f32>()
                .map_err(|_| std::io::ErrorKind::InvalidData)?,
            title_parts[5]
                .parse::<f32>()
                .map_err(|_| std::io::ErrorKind::InvalidData)?,
        );

        Ok(ReadTrack {
            path: full_path,
            time_step: 0,
            file_track,
            num_atoms,
            my_box,
            x: vec![0.0; num_atoms as usize],
            y: vec![0.0; num_atoms as usize],
            z: vec![0.0; num_atoms as usize],
            btype: vec![0; num_atoms as usize],
            iterator: (0..num_atoms as usize).collect(),
        })
    }

    fn one_step(&mut self) -> Result<bool, std::io::Error> {
        let mut title = String::new();
        self.file_track.read_line(&mut title)?;
        let title_parts: Vec<&str> = title.split_whitespace().collect();
        self.time_step = title_parts
            .get(1)
            .ok_or(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "Missing time step",
            ))?
            .parse::<i32>()
            .map_err(|_| {
                std::io::Error::new(std::io::ErrorKind::InvalidData, "Invalid time step")
            })?;
        for i in &self.iterator {
            let mut record = String::new();
            self.file_track.read_line(&mut record)?;
            let record_parts: Vec<&str> = record.split_whitespace().collect();
            self.x[*i] = record_parts[1]
                .parse::<f32>()
                .map_err(|_| std::io::ErrorKind::InvalidData)?;
            self.y[*i] = record_parts[2]
                .parse::<f32>()
                .map_err(|_| std::io::ErrorKind::InvalidData)?;
            self.z[*i] = record_parts[3]
                .parse::<f32>()
                .map_err(|_| std::io::ErrorKind::InvalidData)?;
            self.btype[*i] = record_parts[4]
                .parse::<i32>()
                .map_err(|_| std::io::ErrorKind::InvalidData)?;
        }

        Ok(true)
    }
}

use std::f32::consts::PI;

fn rdf(group_a: &Vec<[f32; 3]>, group_b: &Vec<[f32; 3]>, box_obj: &MyBox, dr: f32) -> Vec<f32> {
    let min_box = box_obj.x.min(box_obj.y).min(box_obj.z);
    let num_bins = (2.0 * min_box / dr).ceil() as usize;
    let mut g = vec![0.0; num_bins];
    let mut v = vec![0.0; num_bins];

    for (j, v_elem) in v.iter_mut().enumerate() {
        *v_elem = (4.0 / 3.0) * PI * (((j + 1) as f32 * dr).powi(3) - (j as f32 * dr).powi(3));
    }

    for atm1 in group_a {
        for atm2 in group_b {
            let dx = atm1[0] - atm2[0];
            let dy = atm1[1] - atm2[1];
            let dz = atm1[2] - atm2[2];
            let (dx, dy, dz) = box_obj.periodic_correct(dx, dy, dz);
            let dist = (dx.powi(2) + dy.powi(2) + dz.powi(2)).sqrt();
            g[(dist / dr).floor() as usize] += 1.0;
        }
    }

    let v_scaled = v
        .iter()
        .map(|&v_elem| v_elem / group_b.len() as f32)
        .collect::<Vec<f32>>();
    g.iter()
        .zip(v_scaled.iter())
        .map(|(&g_elem, &v_elem)| g_elem / v_elem)
        .collect()
}

fn main() -> Result<(), Box<dyn Error>> {
    let root = BitMapBackend::new("rdf.png", (640, 480)).into_drawing_area();
    root.fill(&WHITE)?;
    let mut rt = ReadTrack::new(r"C:\Users\siera\code\new\Scan_Track\brush").unwrap();
    let mut label = true;
    let mut group_a: Vec<[f32; 3]> = Vec::new(); // Change the type to Vec<[f32; 3]>
    let mut group_b: Vec<[f32; 3]> = Vec::new(); // Change the type to Vec<[f32; 3]>
    let mut _data: Vec<(f32, f32)> = Vec::new();
    let mut x_values: Vec<f32> = Vec::new();
    let mut y_values: Vec<f32> = Vec::new();

    while let Ok(result) = rt.one_step() {
        if label {
            let mut n_a = 0;
            let mut n_b = 0;
            for i in &rt.btype {
                if *i == 1 {
                    n_a += 1;
                } else if *i == 2 {
                    n_b += 1;
                }
            }
            group_a = vec![[0.0; 3]; n_a]; // Change the inner type to [f32; 3]
            group_b = vec![[0.0; 3]; n_b]; // Change the inner type to [f32; 3]
            label = false;
        }
        let mut it = 0;
        let mut ik = 0;
        for i in 0..rt.num_atoms as usize {
            if rt.btype[i] == 1 {
                group_a[it] = [rt.x[i], rt.y[i], rt.z[i]]; // Change to array [f32; 3]
                it += 1;
            } else if rt.btype[i] == 2 {
                group_b[ik] = [rt.x[i], rt.y[i], rt.z[i]]; // Change to array [f32; 3]
                ik += 1;
            }
        }
        if !result {
            break; // Break the loop if the result is false
        }
    }
    println!("{:?}", rt.time_step);
    println!("{:?}", group_a);
    println!("{:?}", group_b);

    let dr = 0.1; // Example value for dr, modify as per your needs

    let rdff = rdf(&group_a, &group_b, &rt.my_box, dr); // Provide the missing `dr` argument
    let x_values: Vec<f32> = (0..rdff.len()).map(|i| i as f32 * 0.1).collect();
    let y_values: Vec<f32> = rdff.iter().cloned().collect();

    let root = BitMapBackend::new(
        r"C:\Users\siera\code\new\Scan_Track\brush\rdf.png",
        (800, 600),
    )
    .into_drawing_area();
    root.fill(&WHITE)?;

    let (x_min, x_max) = (x_values[0], x_values[x_values.len() - 1]);
    let (y_min, y_max) = (
        y_values.iter().fold(f32::INFINITY, |a, &b| a.min(b)),
        y_values.iter().fold(f32::NEG_INFINITY, |a, &b| a.max(b)),
    );

    let mut chart = ChartBuilder::on(&root)
        .x_label_area_size(40)
        .y_label_area_size(40)
        .margin(5)
        .caption("RDF Plot", ("sans-serif", 30.0).into_font())
        .build_cartesian_2d(x_min..x_max, y_min..y_max)?;

    chart.configure_mesh().x_labels(10).y_labels(10).draw()?;

    chart.draw_series(LineSeries::new(
        x_values.iter().zip(y_values.iter()).map(|(&x, &y)| (x, y)),
        &RED,
    ))?;
    chart.configure_mesh().x_labels(10).y_labels(10).draw()?;

    chart.draw_series(LineSeries::new(
        x_values.iter().zip(y_values.iter()).map(|(&x, &y)| (x, y)),
        &RED,
    ))?;

    root.present()?;

    let mut file = File::create(r"C:\Users\siera\code\new\Scan_Track\brush\rdf_out.txt")?;
    file.write_all(format!("{:?}", rdff).as_bytes())?;

    Ok(())
}
