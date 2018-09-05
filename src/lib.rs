extern crate image;
extern crate num_traits;
extern crate rustc_serialize;

use std::vec::Vec;
use std::f64;
use std::path::Path;
use std::cmp::{max, min};

use image::{Pixel, GenericImage};

use num_traits::{ToPrimitive};

#[derive(Debug)]
#[derive(Default)]
#[derive(Clone)]
pub struct RGBColor {
    pub r : u8,
    pub g : u8,
    pub b : u8,
}

#[derive(Debug)]
#[derive(Default)]
#[derive(Clone)]
pub struct XYZColor {
    pub x : f64,
    pub y : f64,
    pub z : f64,
}

#[derive(Debug)]
#[derive(Default)]
#[derive(Clone)]
pub struct LABColor {
    pub l : f64,
    pub a : f64,
    pub b : f64,
}

#[derive(Debug)]
#[derive(Default)]
#[derive(Clone)]
pub struct ClusterInfo {
    pub lab : LABColor,
    pub x : f64,
    pub y : f64,
    pub z : f64,
}


fn xyz_pre_transform(val : f64) -> f64
{
    return if val <= 0.04045 {
            val/12.92
        } else {
            ((val+0.055)/1.055).powf(2.4)
        };
}

pub fn rgb_to_xyz (rgb : &RGBColor, xyz : &mut XYZColor)
{
    let RGBColor {r: sr, g: sg, b: sb} = *rgb;
    
    let mut r : f64 = (sr as f64) / 255.0;
    let mut g : f64 = (sg as f64) / 255.0;
    let mut b : f64 = (sb as f64) / 255.0;

    r =	xyz_pre_transform(r);
    g =	xyz_pre_transform(g);
    b =	xyz_pre_transform(b);

    *xyz = XYZColor {
        x: r*0.4124564 + g*0.3575761 + b*0.1804375,
        y: r*0.2126729 + g*0.7151522 + b*0.0721750,
        z: r*0.0193339 + g*0.1191920 + b*0.9503041,
    };
}


pub fn rgb_to_lab(rgb : &RGBColor, lab: &mut LABColor)
{
    let mut xyz : XYZColor = Default::default();
    rgb_to_xyz(rgb, &mut xyz);

    // CIE standard
    let epsilon = 0.008856;
    let kappa   = 903.3;

    // reference white
    let xyz_white = XYZColor{x : 0.950456, y: 1.0, z: 1.088754};
    
    let mut xr = xyz.x / xyz_white.x;
    let mut yr = xyz.y / xyz_white.y;
    let mut zr = xyz.z / xyz_white.z;
    
    xr = if xr > epsilon {
        xr.powf(1.0/3.0)
    } else {
        (kappa*xr + 16.0)/116.0
    };
    yr = if yr > epsilon {
        yr.powf(1.0/3.0)
    } else {
        (kappa*yr + 16.0)/116.0
    };
    zr = if zr > epsilon {
        zr.powf(1.0/3.0)
    } else {
        (kappa*zr + 16.0)/116.0
    };

    *lab = LABColor {
        l: 116.0*yr-16.0,
        a: 500.0*(xr-yr),
        b: 200.0*(yr-zr),
    }
}

pub fn rgb_to_lab_image<I: GenericImage>(
    img: &I,
    lab_vec: &mut Vec<LABColor>
    )
{
    for (x, y, pixel) in (*img).pixels() {
        let sp = pixel.to_rgb();
        let rgb = RGBColor{
            r: sp[0].to_u8().unwrap(),
            g: sp[1].to_u8().unwrap(),
            b: sp[2].to_u8().unwrap()
        };
        let idx = ((y * (*img).width()) + x) as usize;

        rgb_to_lab(&rgb, &mut lab_vec[idx]);
    }
}

fn init_seeds(
    seeds: &mut Vec<ClusterInfo>,
    width: i32, height: i32,
    lab_vec: & Vec<LABColor>,
    step: i32)
{
    
    let s = step as f64;

    let mut xstrips = (0.5 + (width as f64)/s) as i32;
    let mut ystrips = (0.5 + (height as f64)/s) as i32;

    let mut xerr = width - step*xstrips;
    if xerr < 0 {
        xstrips -= 1;
        xerr = width - step*xstrips;
    }

    let mut yerr = height - step*ystrips;
    if yerr < 0 {
        ystrips -= 1;
        yerr = height - step*ystrips;
    }

    let xerrperstrip = (xerr as f64)/(xstrips as f64);
    let yerrperstrip = (yerr as f64)/(ystrips as f64);
    
    let xoff = step/2;
    let yoff = step/2;

    let numseeds = (xstrips*ystrips) as usize;

    for _i in 0..numseeds {
        seeds.push(Default::default());
    }
    
    let mut n = 0usize;
    
    for y in 0..ystrips
    {
        let ye = ((y as f64) * yerrperstrip) as i32;
        
        for x in 0..xstrips
        {
            let xe = ((x as f64) * xerrperstrip) as i32;
            
            let i : usize = ((y*step+yoff+ye)*width + (x*step+xoff+xe)) as usize;
            
            seeds[n] = ClusterInfo {
                lab: lab_vec[i].clone(),
                x: (x*step+xoff+xe) as f64,
                y: (y*step+yoff+ye) as f64,
                z: 0.0,
            };
            n += 1;
        }
    }
}

pub fn slic<I: GenericImage>(
    img: &I,
    lab_vec: & Vec<LABColor>,
    seeds_vec: &mut Vec<ClusterInfo>,
    klabels: &mut Vec<i32>,
    step: i32,
    compactness: f64,
    iterations: u32
    )
{
    let width = img.width() as i32;
    let height = img.height() as i32;
    let depth = 1 as i32;

    let sz : usize = (img.width() * img.height()) as usize;
    
    let numk = seeds_vec.len();
    let offset = step;
    
    let mut cluster_size: Vec<f64> = vec![0.0; numk];
    let mut inv: Vec<f64> = vec![0.0; numk];

    let mut sigma_l: Vec<f64> = vec![0.0; numk];
    let mut sigma_a: Vec<f64> = vec![0.0; numk];
    let mut sigma_b: Vec<f64> = vec![0.0; numk];
    let mut sigma_x: Vec<f64> = vec![0.0; numk];
    let mut sigma_y: Vec<f64> = vec![0.0; numk];
    let mut sigma_z: Vec<f64> = vec![0.0; numk];
    let init_sigma: Vec<f64> = vec![0.0; numk];

    let init_distances: Vec<f64> = vec![f64::MAX; sz];
    let mut dist_vec: Vec<f64> =  vec![f64::MAX; sz];

    let invwt = 1.0 / ((step as f64/compactness)*(step as f64/compactness));

    let mut x1: u32; let mut y1: u32; let mut x2: u32; let mut y2: u32; let mut z1: u32; let mut z2: u32;
    
    let mut dist: f64;
    let mut dist_xyz: f64;
    
    let slic_iterations : i32 = iterations as i32;
    for _i in 0..slic_iterations
    {
        dist_vec.clone_from(&init_distances);
        for n in 0..numk
        {
            let ref s = seeds_vec[n];

            y1 = max(0i32,		s.y as i32-offset) as u32;
            y2 = min(height,	s.y as i32+offset) as u32;
            x1 = max(0i32,		s.x as i32-offset) as u32;
            x2 = min(width,		s.x as i32+offset) as u32;
            z1 = max(0i32,		s.z as i32-offset) as u32;
            z2 = min(depth,		s.z as i32+offset) as u32;

            for z in z1..z2
            {
                for y in y1..y2
                {
                    for x in x1..x2
                    {
                        let i = (y*(width as u32) + x) as usize;

                        
                        let l = lab_vec[i].l;
                        let a = lab_vec[i].a;
                        let b = lab_vec[i].b;
                        
                        let ref c = s.lab;

                        dist =			((l - c.l)*(l - c.l) +
                                        (a - c.a)*(a - c.a) +
                                        (b - c.b)*(b - c.b)) as f64;

                        dist_xyz =		((x as f64 - s.x)*(x as f64 - s.x) +
                                        (y as f64 - s.y)*(y as f64 - s.y) +
                                        (z as f64 - s.z)*(z as f64 - s.z)) as f64;
                        dist += dist_xyz*invwt;

                        //println!("z {} y {} x {} i {} dist {} dist_vec[i] {}", z, y, x, i, dist, dist_vec[i]);

                        if dist < dist_vec[i]
                        {
                            dist_vec[i] = dist;
                            klabels[i] = n as i32;
                        }
                    }
                }
            }
        }
        // Recalculate centroid and store in the seed values
        sigma_l.clone_from(&init_sigma);
        sigma_a.clone_from(&init_sigma);
        sigma_b.clone_from(&init_sigma);
        sigma_x.clone_from(&init_sigma);
        sigma_y.clone_from(&init_sigma);
        sigma_z.clone_from(&init_sigma);
        cluster_size.clone_from(&init_sigma);

        for d in 0..depth
        {
            let mut ind : usize = 0;
            for r in 0..height
            {
                for c in 0..width
                {
                    let l = klabels[ind] as usize;
                    sigma_l[l] += lab_vec[ind].l;
                    sigma_a[l] += lab_vec[ind].a;
                    sigma_b[l] += lab_vec[ind].b;
                    sigma_x[l] += c as f64;
                    sigma_y[l] += r as f64;
                    sigma_z[l] += d as f64;

                    cluster_size[l] += 1.0;
                    ind += 1;
                }
            }
        }

        for k in 0..numk
        {
            if cluster_size[k] <= 0.0 { cluster_size[k] = 1.0 };
            // Computing inverse now to multiply, instead of dividing later
            inv[k] = 1.0/cluster_size[k];
        }
        
        for k in 0..numk
        {
            let mut c = &mut seeds_vec[k];
            c.lab.l = (sigma_l[k]*inv[k]) as f64;
            c.lab.a = (sigma_a[k]*inv[k]) as f64;
            c.lab.b = (sigma_b[k]*inv[k]) as f64;
            c.x = (sigma_x[k]*inv[k]) as f64;
            c.y = (sigma_y[k]*inv[k]) as f64;
            c.z = (sigma_z[k]*inv[k]) as f64;
        }
    }

}

pub fn enforce_label_connectivity(
    labels: &mut Vec<i32>,
    width: i32,
    height: i32,
    nlabels: &mut Vec<i32>,
    numlabels: &mut u32,
    step: &usize
    )
{
    let dx4 : Vec<i32> = vec![-1,  0,  1,  0];
    let dy4 : Vec<i32> = vec![ 0, -1,  0,  1];
    
    //let sz = width*height;
    let supsz = (step * step) as i32;

    let mut adj_label = 0; //adjacent label
    let mut xvec = vec![0; (supsz * 10) as usize];
    let mut yvec = vec![0; (supsz * 10) as usize];
    
    //let mut nlabels = vec![-1; sz as usize];
    
    // Labeling
    let mut lab = 0i32;
    
    let mut i = 0;

    for h in 0..height
    {
        for w in 0..width
        {
            if nlabels[i] < 0
            {
                nlabels[i] = lab;
                // Quickly find an adjacent label for use later if needed
                for n in 0..4
                {
                    let x = w + dx4[n];
                    let y = h + dy4[n];
                    
                    if (x >= 0 && x < width) && (y >= 0 && y < height)
                    {
                        let nindex = (y*width + x) as usize;
                        if nlabels[nindex] >= 0
                        {
                            adj_label = nlabels[nindex];
                        }
                    }
                }
                
                xvec[0] = w;
                yvec[0] = h;
                let mut count = 1;
                let mut c = 0;
                while c < count
                {
                    for n in 0..4
                    {
                        let x = xvec[c] + dx4[n];
                        let y = yvec[c] + dy4[n];
                        
                        if (x >= 0 && x < width) && (y >= 0 && y < height)
                        {
                            let nindex = (y*width + x) as usize;

                            if 0 > nlabels[nindex] && labels[i] == labels[nindex]
                            {
                                xvec[count] = x;
                                yvec[count] = y;
                                nlabels[nindex] = lab;
                                count += 1;
                            }
                        }

                    }
                    c += 1;
                }
                
                // Remove small segments. If segment size is less then a limit, assign
                // adjacent label found previously, then decrement label count.
                if count <= (supsz >> 2) as usize
                {
                    for c in 0..count
                    {
                        let ind = (yvec[c]*width+xvec[c]) as usize;
                        nlabels[ind] = adj_label;
                    }
                    //println!("label {:?} too small at size {:?} (limit {:?})", lab, count, (supsz >> 2) as usize);
                    lab -= 1;
                }
                lab += 1;
            }
            i += 1;
        }
    }
    
    (*numlabels) = lab as u32;
}

pub fn save_image_with_segment_mean_colour(
    img: & image::RgbImage,
    klabels: &Vec<i32>,
    out_fn: &str,
    outline_color: &image::Rgb<u8>
    )
{
    let (width, height) = (img.width() as u32, img.height() as u32);

    // Create a new ImgBuf with width: imgx and height: imgy
    let mut mean_img = image::ImageBuffer::new(width, height);
    

    for x in 0..width {
        for y in 0..height {
            let idx = y * width + x;
            let px = img.get_pixel(x, y);
            mean_img.put_pixel(x, y, *px);
        }
    }
            
    // Write the contents of this image to the Writer in PNG format.
    let _ = mean_img.save(&Path::new(out_fn)).unwrap();

}

pub fn save_image_with_segment_boundaries(
    img: & image::RgbImage,
    klabels: &Vec<i32>,
    out_fn: &str,
    outline_color: &image::Rgb<u8>
    )
{
    let (width, height) = (img.width() as u32, img.height() as u32);

    // Create a new ImgBuf with width: imgx and height: imgy
    let mut contour = image::ImageBuffer::new(width, height);
    

    for x in 0..width {
        for y in 0..height {
            let idx = y * width + x;
            // only if y > 0
            let up_idx = if y > 0 {
                ((y - 1) * width + x)
            } else {
                idx
            };

            // only if x > 0
            let left_idx = if x > 0 {
                (y * width + (x - 1))
            } else {
                idx
            };

            if klabels[up_idx as usize] != klabels[idx as usize] {
                contour.put_pixel(x, y, *outline_color);
            } else if klabels[left_idx as usize] != klabels[idx as usize] {
                contour.put_pixel(x, y, *outline_color);
            } else {
                let px = img.get_pixel(x, y);
                contour.put_pixel(x, y, *px);
            }
        }
    }
            
    // Write the contents of this image to the Writer in PNG format.
    let _ = contour.save(&Path::new(out_fn)).unwrap();
}

pub fn do_segmentation_with_superpixel_size(
    img: & image::RgbImage,
    lab_vec: &mut Vec<LABColor>,
    width: u32,
    height: u32,
    klabels: &mut Vec<i32>,
    numlabels: &mut u32,
    superpixel_size: u32, 
    compactness: f64,
    iterations: u32 
    ) -> u32 // returns actual number of superpixels
{
    //------------------------------------------------
    let step: u32 = ((superpixel_size as f64).sqrt()+0.5) as u32;
    
    let mut seeds_vec : Vec<ClusterInfo> = vec![];

    let sz : usize = (width*height) as usize;
    
    klabels.clear();
    for _s in 0..sz {
        klabels.push(-1);	
    } 

    rgb_to_lab_image(img, lab_vec);
    
    init_seeds(&mut seeds_vec, img.width() as i32, img.height() as i32, &lab_vec, step as i32);

    slic(img, lab_vec, &mut seeds_vec, klabels, step as i32, compactness, iterations);
    
    (*numlabels) = seeds_vec.len() as u32;

    let mut new_labels : Vec<i32> = vec![];
    for _s in 0..sz {
        new_labels.push(-1);	
    } 

    enforce_label_connectivity(klabels, width as i32, height as i32, &mut new_labels, numlabels, &(step as usize));
    for x in 0usize..(sz as usize) {
        klabels[x] = new_labels[x];
    }
    //let outline_color = image::Rgb([255 as u8, 255 as u8, 255 as u8]);

    return *numlabels;
}

pub fn do_segmentation_with_num_superpixels(
    img: & image::RgbImage,
    lab_vec: &mut Vec<LABColor>,
    width: u32,
    height: u32,
    klabels: &mut Vec<i32>,
    num_superpixels: u32, 
    compactness: f64,
    iterations: u32 
    ) -> u32 // returns actual number of superpixels
{
    let mut numlabels: u32 = 0;
    let superpixel_size : f64 = 0.5 + (width*height) as f64 / num_superpixels as f64;
    return do_segmentation_with_superpixel_size(
        img, lab_vec, width, height, klabels, &mut numlabels, superpixel_size as u32, compactness, iterations);
}

#[cfg(test)]
mod test 
{
    use std::f64;

    fn almost_eq(a : f64, b : f64)
    {
        let diff = (a - b).abs();
        assert!(diff < 0.0001, format!("{} and {} are not within 0.0001", a, b));
    }

    #[test]
    fn test_rgb_to_xyz()
    {
        let rgb = super::RGBColor{r: 255, g: 255, b: 22};
        let mut xyz : super::XYZColor = Default::default();
        super::rgb_to_xyz(&rgb, &mut xyz);

        almost_eq(xyz.x, 0.7714801);
        almost_eq(xyz.y, 0.9284041);
        almost_eq(xyz.z, 0.1461503);
    }

    #[test]
    fn test_rgb_to_lab()
    {
        let rgb = super::RGBColor{r: 255, g: 255, b: 22};
        let mut lab : super::LABColor = Default::default();
        super::rgb_to_lab(&rgb, &mut lab);

        almost_eq(lab.l, 97.1627998);
        almost_eq(lab.a, -21.3609707);
        almost_eq(lab.b, 92.7035375);
    }
}
