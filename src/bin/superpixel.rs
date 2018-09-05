extern crate superpixel;
extern crate image;
extern crate num_traits;
extern crate rustc_serialize;
extern crate docopt;

use std::vec::Vec;
use std::path::Path;
use std::iter::repeat;
use std::time::Instant;

use image::DynamicImage;

use superpixel::*;

use docopt::Docopt;

const USAGE: &'static str = "
superpixel

Usage:
  superpixel [--compactness=<compactness>] [--patches=<patches>] [--iterations=<iterations>] [--outline-image=<outline_image> ] <image>
  superpixel (-h | --help)
  superpixel --version

Options:
  -h --help       Show this screen.
  --version       Show version.
  --compactness=<compactness>     Specify compactness parameter. [default: 20.0]
  --patches=<patches>     Specify number of patches to split image into. [default: 8000]
  --iterations=<iterations>     Number of iterations to run k-means. [default: 10]
  --outline-image=<outline_image>  Image file to write segment visualisation to. [default: segment_vis.png]
  <image>   	  Image file to segment.
";

#[derive(Debug, RustcDecodable)]
struct Args {
    pub flag_compactness: Option<f64>,
    pub flag_patches: Option<u32>,
    pub flag_iterations: Option<u32>,
    pub flag_outline_image: Option<String>,
    pub arg_image: String,
}

fn main() {
    let args: Args = Docopt::new(USAGE)
                            .and_then(|d| d.decode())
                            .unwrap_or_else(|e| e.exit());


	// Use the open function to load an image from a Path.
    // ```open``` returns a dynamic image.
    let x = image::open(&Path::new(&args.arg_image[..])).unwrap();
    let img = match x {
    	DynamicImage::ImageRgb8(val) =>
    	 	val,
    	_ =>
    		unreachable!()
    	};

    // The dimensions method returns the images width and height
    println!("Dimensions: width {:?}px, height {:?}px", img.dimensions().0, img.dimensions().1);

    let mut lab_vec : Vec<LABColor> = repeat(Default::default()).take((img.width() * img.height()) as usize).collect();
    let now = Instant::now();
    print!("Converting to LAB color space... ");
    rgb_to_lab_image(&img, &mut lab_vec);
    let new_now = Instant::now();
    let elapsed = new_now.duration_since(now);
    println!("{:?}.{:?}s", elapsed.as_secs(), elapsed.subsec_millis());

    let mut klabels : Vec<i32> = vec![];
    let now = Instant::now();
    print!("Segmenting... ");
    let num_superpixels = do_segmentation_with_num_superpixels(
        &img, &mut lab_vec,
        img.width(), img.height(),
        &mut klabels,
        args.flag_patches.unwrap(), args.flag_compactness.unwrap(), args.flag_iterations.unwrap());
    let new_now = Instant::now();
    let elapsed = new_now.duration_since(now);
    print!("{:?}.{:?}s ", elapsed.as_secs(), elapsed.subsec_millis());
    println!("(with {} segments)", num_superpixels);

    let outline_fn = args.flag_outline_image.unwrap();
    println!("Saving outline image to {:?}", outline_fn);
    let outline_color = image::Rgb([255u8, 255u8, 255u8]);
    save_image_with_segment_boundaries(&img, &klabels, &outline_fn, &outline_color);
}

