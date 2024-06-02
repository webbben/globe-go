//! Customizable ASCII globe generator.
//!
//! Based on [C++ code by DinoZ1729](https://github.com/DinoZ1729/Earth).

#![allow(dead_code)]

use std::f32::consts::PI;
use std::fs::File;

use std::io::{stdout, Read, Stdout, Write};
use std::time::Duration;

use crossterm::{
    cursor,
    event::{poll, read, Event, KeyCode},
    style::Print,
    ExecutableCommand, QueueableCommand,
};
use crossterm::{event::MouseEvent, terminal};

use crossterm::terminal::ClearType;

#[no_mangle]
pub extern "C" fn ext_screensaver(settings_ffi: SettingsFFI) {
    // convert to regular Settings struct
    let settings = Settings {
        refresh_rate: settings_ffi.refresh_rate,
        globe_rotation_speed: settings_ffi.globe_rotation_speed,
        cam_rotation_speed: settings_ffi.cam_rotation_speed,
        cam_zoom: settings_ffi.cam_zoom,
        focus_speed: settings_ffi.focus_speed,
        night: settings_ffi.night,
        coords: (settings_ffi.coord_x, settings_ffi.coord_y),
    };

    //start_screensaver(settings);
    my_screensaver(settings)
}

#[no_mangle]
pub extern "C" fn ext_interactive(settings_ffi: SettingsFFI) {
    // convert to regular Settings struct
    let settings = Settings {
        refresh_rate: settings_ffi.refresh_rate,
        globe_rotation_speed: settings_ffi.globe_rotation_speed,
        cam_rotation_speed: settings_ffi.cam_rotation_speed,
        cam_zoom: settings_ffi.cam_zoom,
        focus_speed: settings_ffi.focus_speed,
        night: settings_ffi.night,
        coords: (settings_ffi.coord_x, settings_ffi.coord_y),
    };
    start_interactive(settings)
}

fn my_screensaver(settings: Settings) {
    let mut term_size = terminal::size().unwrap();

    let mut globe = GlobeConfig::new()
        //.use_template(GlobeTemplate::Earth)
        .with_texture(EARTH_NIGHT_TEXTURE, None)
        .with_camera(CameraConfig::default())
        .display_night(settings.night)
        .build();

    // create a new canvas
    //let mut canvas = Canvas::new(250, 250, None);
    let mut canvas = if term_size.0 > term_size.1 {
        Canvas::new(term_size.1 * 8, term_size.1 * 8, None)
    } else {
        Canvas::new(term_size.0 * 4, term_size.0 * 4, None)
    };
    let mut cam_xy = 0.;
    let mut cam_z = 0.;

    loop {
        if poll(Duration::from_millis(1000 / settings.refresh_rate as u64)).unwrap() {
            match read().unwrap() {
                // pressing any key exists the program
                Event::Key(_) => break,
                Event::Resize(width, height) => {
                    term_size = (width, height);
                    canvas = if width > height {
                        Canvas::new(height * 8, height * 8, None)
                    } else {
                        Canvas::new(width * 4, width * 4, None)
                    };
                }
                Event::Mouse(_) => (),
            }
        }

        // apply globe rotation
        globe.angle += settings.globe_rotation_speed;
        //cam_xy -= settings.globe_rotation_speed / 2.;

        // apply camera rotation
        //cam_xy -= settings.cam_rotation_speed;

        //globe.camera.update(settings.cam_zoom, cam_xy, cam_z);

        // render globe on the canvas
        canvas.clear();
        globe.render_on(&mut canvas);

        // print out the canvas
        let (size_x, size_y) = canvas.get_size();
        for i in 0..size_y / 8 {
            for j in 0..size_x / 4 {
                print!("{}", canvas.matrix[i][j]);
            }
            println!();
        }
    }
}

pub type Int = isize;
pub type Float = f32;

static EARTH_TEXTURE: &str = include_str!("../textures/earth.txt");
static EARTH_NIGHT_TEXTURE: &str = include_str!("../textures/earth_night.txt");

/// Collection of scene settings that get passed from clap to mode processing
/// functions.
pub struct Settings {
    /// Refresh rate in cycles per second
    refresh_rate: usize,
    /// Initial globe rotation speed
    globe_rotation_speed: f32,
    /// Initial camera rotation speed
    cam_rotation_speed: f32,
    /// Initial camera zoom
    cam_zoom: f32,
    /// Target focus speed
    focus_speed: f32,
    /// Globe night side switch
    night: bool,
    /// Initial location coordinates
    coords: (f32, f32),
}

// Settings struct that is FFI compatible
#[repr(C)]
pub struct SettingsFFI {
    /// Refresh rate in cycles per second
    refresh_rate: usize,
    /// Initial globe rotation speed
    globe_rotation_speed: f32,
    /// Initial camera rotation speed
    cam_rotation_speed: f32,
    /// Initial camera zoom
    cam_zoom: f32,
    /// Target focus speed
    focus_speed: f32,
    /// Globe night side switch
    night: bool,
    /// Initial location coordinates
    coord_x: f32,
    coord_y: f32,
}

/// Globe texture.
pub struct Texture {
    day: Vec<Vec<char>>,
    night: Option<Vec<Vec<char>>>,
    palette: Option<Vec<char>>,
}

impl Texture {
    pub fn new(
        day: Vec<Vec<char>>,
        night: Option<Vec<Vec<char>>>,
        palette: Option<Vec<char>>,
    ) -> Self {
        Texture {
            day,
            night,
            palette,
        }
    }
    pub fn get_size(&self) -> (usize, usize) {
        (self.day[0].len() - 1, self.day.len() - 1)
    }
}

/// Canvas that will be used to render the globe onto.
pub struct Canvas {
    pub matrix: Vec<Vec<char>>,
    size: (usize, usize),
    // character size
    pub char_pix: (usize, usize),
}

impl Canvas {
    pub fn new(x: u16, y: u16, cp: Option<(usize, usize)>) -> Self {
        let x = x as usize;
        let y = y as usize;

        let matrix = vec![vec![' '; x]; y];

        Self {
            size: (x, y),
            matrix,
            char_pix: cp.unwrap_or((4, 8)),
        }
    }
    pub fn get_size(&self) -> (usize, usize) {
        self.size
    }
    pub fn clear(&mut self) {
        for i in self.matrix.iter_mut().flatten() {
            *i = ' ';
        }
    }
    fn draw_point(&mut self, a: usize, b: usize, c: char) {
        if a >= self.size.0 || b >= self.size.1 {
            return;
        }
        self.matrix[b][a] = c;
    }
}

/// Main globe abstraction.
pub struct Globe {
    pub camera: Camera,
    pub radius: Float,
    pub angle: Float,
    pub texture: Texture,
    pub display_night: bool,
}

impl Globe {
    pub fn render_on(&self, canvas: &mut Canvas) {
        // let there be light
        let light: [Float; 3] = [0., 999999., 0.];
        // shoot the ray through every pixel
        let (size_x, size_y) = canvas.get_size();
        for yi in 0..size_y {
            let yif = yi as Int;
            for xi in 0..size_x {
                let xif = xi as Int;
                // coordinates of the camera, origin of the ray
                let o: [Float; 3] = [self.camera.x, self.camera.y, self.camera.z];
                // u is unit vector, direction of the ray
                let mut u: [Float; 3] = [
                    -((xif - (size_x / canvas.char_pix.0 / 2) as Int) as Float + 0.5)
                        / (size_x / canvas.char_pix.0 / 2) as Float,
                    ((yif - (size_y / canvas.char_pix.1 / 2) as Int) as Float + 0.5)
                        / (size_y / canvas.char_pix.1 / 2) as Float,
                    -1.,
                ];
                transform_vector(&mut u, self.camera.matrix);
                u[0] -= self.camera.x;
                u[1] -= self.camera.y;
                u[2] -= self.camera.z;
                normalize(&mut u);
                let dot_uo = dot(&u, &o);
                let discriminant: Float = dot_uo * dot_uo - dot(&o, &o) + self.radius * self.radius;

                // ray doesn't hit the sphere
                if discriminant < 0. {
                    continue;
                }

                let distance: Float = -discriminant.sqrt() - dot_uo;

                // intersection point
                let inter: [Float; 3] = [
                    o[0] + distance * u[0],
                    o[1] + distance * u[1],
                    o[2] + distance * u[2],
                ];

                // surface normal
                let mut n: [Float; 3] = [
                    o[0] + distance * u[0],
                    o[1] + distance * u[1],
                    o[2] + distance * u[2],
                ];
                normalize(&mut n);

                // unit vector pointing from intersection to light source
                let mut l: [Float; 3] = [0.; 3];
                vector(&mut l, &inter, &light);
                normalize(&mut l);
                let luminance: Float = clamp(5. * (dot(&n, &l)) + 0.5, 0., 1.);
                let mut temp: [Float; 3] = [inter[0], inter[1], inter[2]];
                rotate_x(&mut temp, -PI * 2. * 0. / 360.);

                // computing coordinates for the sphere
                let phi: Float = -temp[2] / self.radius / 2. + 0.5;
                let mut theta: Float = (temp[1] / temp[0]).atan() / PI + 0.5 + self.angle / 2. / PI;
                // let mut theta: Float = (temp[1] / temp[0]).atan() / PI + self.angle / 2. / PI * 20.;
                theta -= theta.floor();
                let (tex_x, tex_y) = self.texture.get_size();
                let earth_x = (theta * tex_x as Float) as usize;
                let earth_y = (phi * tex_y as Float) as usize;

                // if night texture and palette are available, draw the night side
                if self.display_night
                    && self.texture.night.is_some()
                    && self.texture.palette.is_some()
                {
                    let palette = self.texture.palette.as_ref().unwrap();
                    let day = find_index(self.texture.day[earth_y][earth_x], palette);
                    let night = find_index(
                        self.texture.night.as_ref().unwrap()[earth_y][earth_x],
                        palette,
                    );

                    let mut index =
                        ((1.0 - luminance) * night as Float + luminance * day as Float) as usize;
                    if index >= palette.len() {
                        index = 0;
                    }
                    canvas.draw_point(xi, yi, palette[index]);
                }
                // else just draw the day texture without considering luminance
                else {
                    canvas.draw_point(xi, yi, self.texture.day[earth_y][earth_x]);
                }
            }
        }
    }
}

/// Globe configuration struct implementing the builder pattern.
#[derive(Default)]
pub struct GlobeConfig {
    camera_cfg: Option<CameraConfig>,
    radius: Option<Float>,
    angle: Option<Float>,
    template: Option<GlobeTemplate>,
    texture: Option<Texture>,
    display_night: bool,
}

impl GlobeConfig {
    /// Creates an empty `GlobeConfig`.
    pub fn new() -> Self {
        Self::default()
    }

    /// Sets `CameraConfig` to be used by the builder.
    pub fn with_camera(mut self, config: CameraConfig) -> Self {
        self.camera_cfg = Some(config);
        self
    }

    /// Sets the globe radius.
    pub fn with_radius(mut self, r: Float) -> Self {
        self.radius = Some(r);
        self
    }

    /// Selects a template to be used by the builder.
    pub fn use_template(mut self, t: GlobeTemplate) -> Self {
        self.template = Some(t);
        self
    }

    /// Sets the day texture to be displayed on the globe.
    pub fn with_texture(mut self, texture: &str, palette: Option<Vec<char>>) -> Self {
        let mut day = Vec::new();
        let lines = texture.lines();
        for line in lines {
            let row: Vec<char> = line.chars().rev().collect();
            day.push(row);
        }
        if let Some(texture) = &mut self.texture {
            texture.day = day;
        } else {
            self.texture = Some(Texture::new(day, None, palette));
        }
        self
    }

    /// Sets the night texture to be displayed on the globe.
    pub fn with_night_texture(mut self, texture: &str, palette: Option<Vec<char>>) -> Self {
        let mut night = Vec::new();
        let lines = texture.lines();
        for line in lines {
            let row: Vec<char> = line.chars().rev().collect();
            night.push(row);
        }

        if let Some(texture) = &mut self.texture {
            texture.night = Some(night);
        } else {
            self.texture = Some(Texture::new(night.clone(), Some(night), palette));
        }

        self
    }

    /// Sets the day texture to be loaded from the given path.
    pub fn with_texture_at(self, path: &str, palette: Option<Vec<char>>) -> Self {
        let mut file = File::open(path).unwrap();
        let mut out_string = String::new();
        file.read_to_string(&mut out_string).unwrap();
        self.with_texture(&out_string, palette)
    }

    /// Sets the night display toggle to the given value.
    pub fn display_night(mut self, b: bool) -> Self {
        self.display_night = b;
        self
    }

    /// Builds new `Globe` from the collected configuration settings.
    pub fn build(mut self) -> Globe {
        if let Some(template) = &self.template {
            match template {
                GlobeTemplate::Earth => {
                    let palette = vec![
                        ' ', '.', ':', ';', '\'', ',', 'w', 'i', 'o', 'g', 'O', 'L', 'X', 'H', 'W',
                        'Y', 'V', '@',
                    ];
                    self = self
                        .with_texture(EARTH_TEXTURE, Some(palette.clone()))
                        .with_night_texture(EARTH_NIGHT_TEXTURE, Some(palette))
                }
            }
        }
        let texture = self.texture.expect("texture not provided");
        let camera = self
            .camera_cfg
            .unwrap_or_else(CameraConfig::default)
            .build();
        Globe {
            camera,
            radius: self.radius.unwrap_or(1.),
            angle: self.angle.unwrap_or(0.),
            texture,
            display_night: self.display_night,
        }
    }
}

/// Built-in globe template enumeration.
pub enum GlobeTemplate {
    Earth,
    // Moon,
    // Mars,
}

/// Camera configuration struct implementing the builder pattern.
pub struct CameraConfig {
    radius: Float,
    alpha: Float,
    beta: Float,
}

impl CameraConfig {
    /// Creates a new `CameraConfig`.
    ///
    /// # Arguments
    ///
    /// - `r` is the distance from the camera to the origin.
    /// - `alfa` is camera's angle along the xy plane.
    /// - `beta` is camera's angle along z axis.
    pub fn new(radius: Float, alpha: Float, beta: Float) -> Self {
        Self {
            radius,
            alpha,
            beta,
        }
    }

    /// Creates a new `CameraConfig` using default values.
    pub fn default() -> Self {
        Self {
            radius: 2.,
            alpha: 0.,
            beta: 0.,
        }
    }

    /// Builds a camera from the collected config information.
    pub fn build(&self) -> Camera {
        let mut camera = Camera::default();
        camera.update(self.radius, self.alpha, self.beta);
        camera
    }
}

#[derive(Default)]
pub struct Camera {
    x: Float,
    y: Float,
    z: Float,
    matrix: [Float; 16],
    inv: [Float; 16],
}

impl Camera {
    /// Updates the camera using new data.
    pub fn update(&mut self, r: Float, alpha: Float, beta: Float) {
        let sin_a = alpha.sin();
        let cos_a = alpha.cos();
        let sin_b = beta.sin();
        let cos_b = beta.cos();

        let x = r * cos_a * cos_b;
        let y = r * sin_a * cos_b;
        let z = r * sin_b;

        let mut matrix = [0.; 16];

        // matrix
        matrix[3] = 0.;
        matrix[7] = 0.;
        matrix[11] = 0.;
        matrix[15] = 1.;
        // x
        matrix[0] = -sin_a;
        matrix[1] = cos_a;
        matrix[2] = 0.;
        // y
        matrix[4] = cos_a * sin_b;
        matrix[5] = sin_a * sin_b;
        matrix[6] = -cos_b;
        // z
        matrix[8] = cos_a * cos_b;
        matrix[9] = sin_a * cos_b;
        matrix[10] = sin_b;

        matrix[12] = x;
        matrix[13] = y;
        matrix[14] = z;

        let mut inv = [0.; 16];
        invert(&mut inv, matrix);

        self.x = x;
        self.y = y;
        self.z = z;
        self.matrix = matrix;
        self.inv = inv;
    }
}

/// Get index of the given character on the palette.
fn find_index(target: char, palette: &[char]) -> Int {
    for (i, &ch) in palette.iter().enumerate() {
        if target == ch {
            return i as Int;
        }
    }
    -1
}

fn transform_vector(vec: &mut [Float; 3], m: [Float; 16]) {
    let tx: Float = vec[0] * m[0] + vec[1] * m[4] + vec[2] * m[8] + m[12];
    let ty: Float = vec[0] * m[1] + vec[1] * m[5] + vec[2] * m[9] + m[13];
    let tz: Float = vec[0] * m[2] + vec[1] * m[6] + vec[2] * m[10] + m[14];
    vec[0] = tx;
    vec[1] = ty;
    vec[2] = tz;
}

fn invert(inv: &mut [Float; 16], matrix: [Float; 16]) {
    inv[0] = matrix[5] * matrix[10] * matrix[15]
        - matrix[5] * matrix[11] * matrix[14]
        - matrix[9] * matrix[6] * matrix[15]
        + matrix[9] * matrix[7] * matrix[14]
        + matrix[13] * matrix[6] * matrix[11]
        - matrix[13] * matrix[7] * matrix[10];

    inv[4] = -matrix[4] * matrix[10] * matrix[15]
        + matrix[4] * matrix[11] * matrix[14]
        + matrix[8] * matrix[6] * matrix[15]
        - matrix[8] * matrix[7] * matrix[14]
        - matrix[12] * matrix[6] * matrix[11]
        + matrix[12] * matrix[7] * matrix[10];

    inv[8] = matrix[4] * matrix[9] * matrix[15]
        - matrix[4] * matrix[11] * matrix[13]
        - matrix[8] * matrix[5] * matrix[15]
        + matrix[8] * matrix[7] * matrix[13]
        + matrix[12] * matrix[5] * matrix[11]
        - matrix[12] * matrix[7] * matrix[9];

    inv[12] = -matrix[4] * matrix[9] * matrix[14]
        + matrix[4] * matrix[10] * matrix[13]
        + matrix[8] * matrix[5] * matrix[14]
        - matrix[8] * matrix[6] * matrix[13]
        - matrix[12] * matrix[5] * matrix[10]
        + matrix[12] * matrix[6] * matrix[9];

    inv[1] = -matrix[1] * matrix[10] * matrix[15]
        + matrix[1] * matrix[11] * matrix[14]
        + matrix[9] * matrix[2] * matrix[15]
        - matrix[9] * matrix[3] * matrix[14]
        - matrix[13] * matrix[2] * matrix[11]
        + matrix[13] * matrix[3] * matrix[10];

    inv[5] = matrix[0] * matrix[10] * matrix[15]
        - matrix[0] * matrix[11] * matrix[14]
        - matrix[8] * matrix[2] * matrix[15]
        + matrix[8] * matrix[3] * matrix[14]
        + matrix[12] * matrix[2] * matrix[11]
        - matrix[12] * matrix[3] * matrix[10];

    inv[9] = -matrix[0] * matrix[9] * matrix[15]
        + matrix[0] * matrix[11] * matrix[13]
        + matrix[8] * matrix[1] * matrix[15]
        - matrix[8] * matrix[3] * matrix[13]
        - matrix[12] * matrix[1] * matrix[11]
        + matrix[12] * matrix[3] * matrix[9];

    inv[13] = matrix[0] * matrix[9] * matrix[14]
        - matrix[0] * matrix[10] * matrix[13]
        - matrix[8] * matrix[1] * matrix[14]
        + matrix[8] * matrix[2] * matrix[13]
        + matrix[12] * matrix[1] * matrix[10]
        - matrix[12] * matrix[2] * matrix[9];

    inv[2] = matrix[1] * matrix[6] * matrix[15]
        - matrix[1] * matrix[7] * matrix[14]
        - matrix[5] * matrix[2] * matrix[15]
        + matrix[5] * matrix[3] * matrix[14]
        + matrix[13] * matrix[2] * matrix[7]
        - matrix[13] * matrix[3] * matrix[6];

    inv[6] = -matrix[0] * matrix[6] * matrix[15]
        + matrix[0] * matrix[7] * matrix[14]
        + matrix[4] * matrix[2] * matrix[15]
        - matrix[4] * matrix[3] * matrix[14]
        - matrix[12] * matrix[2] * matrix[7]
        + matrix[12] * matrix[3] * matrix[6];

    inv[10] = matrix[0] * matrix[5] * matrix[15]
        - matrix[0] * matrix[7] * matrix[13]
        - matrix[4] * matrix[1] * matrix[15]
        + matrix[4] * matrix[3] * matrix[13]
        + matrix[12] * matrix[1] * matrix[7]
        - matrix[12] * matrix[3] * matrix[5];

    inv[14] = -matrix[0] * matrix[5] * matrix[14]
        + matrix[0] * matrix[6] * matrix[13]
        + matrix[4] * matrix[1] * matrix[14]
        - matrix[4] * matrix[2] * matrix[13]
        - matrix[12] * matrix[1] * matrix[6]
        + matrix[12] * matrix[2] * matrix[5];

    inv[3] = -matrix[1] * matrix[6] * matrix[11]
        + matrix[1] * matrix[7] * matrix[10]
        + matrix[5] * matrix[2] * matrix[11]
        - matrix[5] * matrix[3] * matrix[10]
        - matrix[9] * matrix[2] * matrix[7]
        + matrix[9] * matrix[3] * matrix[6];

    inv[7] = matrix[0] * matrix[6] * matrix[11]
        - matrix[0] * matrix[7] * matrix[10]
        - matrix[4] * matrix[2] * matrix[11]
        + matrix[4] * matrix[3] * matrix[10]
        + matrix[8] * matrix[2] * matrix[7]
        - matrix[8] * matrix[3] * matrix[6];

    inv[11] = -matrix[0] * matrix[5] * matrix[11]
        + matrix[0] * matrix[7] * matrix[9]
        + matrix[4] * matrix[1] * matrix[11]
        - matrix[4] * matrix[3] * matrix[9]
        - matrix[8] * matrix[1] * matrix[7]
        + matrix[8] * matrix[3] * matrix[5];

    inv[15] = matrix[0] * matrix[5] * matrix[10]
        - matrix[0] * matrix[6] * matrix[9]
        - matrix[4] * matrix[1] * matrix[10]
        + matrix[4] * matrix[2] * matrix[9]
        + matrix[8] * matrix[1] * matrix[6]
        - matrix[8] * matrix[2] * matrix[5];

    let mut det: Float =
        matrix[0] * inv[0] + matrix[1] * inv[4] + matrix[2] * inv[8] + matrix[3] * inv[12];

    det = 1.0 / det;

    for inv_i in inv.iter_mut() {
        *inv_i *= det;
    }
}

fn cross(r: &mut [Float; 3], a: [Float; 3], b: [Float; 3]) {
    r[0] = a[1] * b[2] - a[2] * b[1];
    r[1] = a[2] * b[0] - a[0] * b[2];
    r[2] = a[0] * b[1] - a[1] * b[0];
}

fn magnitude(r: &[Float; 3]) -> Float {
    dot(r, r).sqrt()
}

fn normalize(r: &mut [Float; 3]) {
    let len: Float = magnitude(r);
    r[0] /= len;
    r[1] /= len;
    r[2] /= len;
}

fn dot(a: &[Float; 3], b: &[Float; 3]) -> Float {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn vector(a: &mut [Float; 3], b: &[Float; 3], c: &[Float; 3]) {
    a[0] = b[0] - c[0];
    a[1] = b[1] - c[1];
    a[2] = b[2] - c[2];
}

fn transform_vector2(vec: &mut [Float; 3], m: &[Float; 9]) {
    vec[0] = m[0] * vec[0] + m[1] * vec[1] + m[2] * vec[2];
    vec[1] = m[3] * vec[0] + m[4] * vec[1] + m[5] * vec[2];
    vec[2] = m[6] * vec[0] + m[7] * vec[1] + m[8] * vec[2];
}

fn rotate_x(vec: &mut [Float; 3], theta: Float) {
    let a = theta.sin();
    let b = theta.cos();
    let m: [Float; 9] = [1., 0., 0., 0., b, -a, 0., a, b];
    transform_vector2(vec, &m);
}

fn rotate_y(vec: &mut [Float; 3], theta: Float) {
    let a = theta.sin();
    let b = theta.cos();
    let m: [Float; 9] = [b, 0., a, 0., 1., 0., -a, 0., b];
    transform_vector2(vec, &m);
}

fn rotate_z(vec: &mut [Float; 3], theta: Float) {
    let a = theta.sin();
    let b = theta.cos();
    let m: [Float; 9] = [b, -a, 0., a, b, 0., 0., 0., 1.];
    transform_vector2(vec, &m);
}

fn clamp(mut x: Float, min: Float, max: Float) -> Float {
    if x < min {
        x = min;
    } else if x > max {
        x = max;
    }
    x
}

/// Screensaver mode doesn't allow for user input. Any key press exits the
/// program.
fn start_screensaver(settings: Settings) {
    terminal::enable_raw_mode().unwrap();
    let mut stdout = stdout();
    stdout.execute(cursor::Hide).unwrap();
    stdout.execute(cursor::DisableBlinking).unwrap();

    let mut term_size = terminal::size().unwrap();

    let mut canvas = if term_size.0 > term_size.1 {
        Canvas::new(term_size.1 * 8, term_size.1 * 8, None)
    } else {
        Canvas::new(term_size.0 * 4, term_size.0 * 4, None)
    };

    let cam_zoom = settings.cam_zoom;
    let mut cam_xy = 0.;
    let mut cam_z = 0.;

    // set the initial coordinates
    //focus_target(settings.coords, 0., &mut cam_xy, &mut cam_z);

    let mut globe = GlobeConfig::new()
        .use_template(GlobeTemplate::Earth)
        .with_camera(CameraConfig::new(cam_zoom, cam_xy, cam_z))
        .display_night(settings.night)
        .build();

    let globe_rot_speed = settings.globe_rotation_speed / 1000.;
    let cam_rot_speed = settings.cam_rotation_speed / 1000.;

    loop {
        if poll(Duration::from_millis(1000 / settings.refresh_rate as u64)).unwrap() {
            match read().unwrap() {
                // pressing any key exists the program
                Event::Key(_) => break,
                Event::Resize(width, height) => {
                    term_size = (width, height);
                    canvas = if width > height {
                        Canvas::new(height * 8, height * 8, None)
                    } else {
                        Canvas::new(width * 4, width * 4, None)
                    };
                }
                Event::Mouse(_) => (),
            }
        }

        // apply globe rotation
        globe.angle += globe_rot_speed;
        cam_xy -= globe_rot_speed / 2.;

        // apply camera rotation
        cam_xy -= cam_rot_speed;

        //globe.camera.update(cam_zoom, cam_xy, cam_z);

        // render globe on the canvas
        canvas.clear();
        globe.render_on(&mut canvas);

        // print canvas to terminal
        print_canvas(&mut canvas, &term_size, &mut stdout);
    }

    stdout.execute(cursor::Show).unwrap();
    stdout.execute(cursor::EnableBlinking).unwrap();

    terminal::disable_raw_mode().unwrap();
    stdout.execute(terminal::Clear(ClearType::All)).unwrap();
}

/// Orients the camera so that it focuses on the given target coordinates.
pub fn focus_target(coords: (f32, f32), xy_offset: f32, cam_xy: &mut f32, cam_z: &mut f32) {
    let (cx, cy) = coords;
    *cam_xy = (cx * PI) * -1. - 1.5 - xy_offset;
    *cam_z = cy * 3. - 1.5;
}

/// Prints globe canvas to stdout.
fn print_canvas(canvas: &mut Canvas, term_size: &(u16, u16), stdout: &mut Stdout) {
    let (canvas_size_x, canvas_size_y) = canvas.get_size();
    for i in 0..canvas_size_y / canvas.char_pix.1 {
        stdout
            .queue(terminal::Clear(terminal::ClearType::CurrentLine))
            .unwrap();
        for j in 0..canvas_size_x / canvas.char_pix.0 {
            stdout.queue(Print(canvas.matrix[i][j])).unwrap();
        }
        stdout.queue(cursor::MoveDown(1)).unwrap();
        stdout
            .queue(cursor::MoveLeft((canvas_size_x / 4) as u16))
            .unwrap();
        stdout.flush().unwrap();
    }

    if term_size.0 / 2 > term_size.1 {
        stdout
            .execute(crossterm::cursor::MoveTo(
                (canvas_size_x / canvas.char_pix.1) as u16
                    - ((canvas_size_x / canvas.char_pix.1) / canvas.char_pix.0) as u16,
                0,
            ))
            .unwrap();
    }
}

/// Interactive mode allows using mouse and/or keyboard to control the globe.
fn start_interactive(settings: Settings) {
    terminal::enable_raw_mode().unwrap();
    let mut stdout = stdout();
    stdout.execute(cursor::Hide).unwrap();
    stdout.execute(cursor::DisableBlinking).unwrap();
    stdout
        .execute(crossterm::event::EnableMouseCapture)
        .unwrap();

    let mut term_size = terminal::size().unwrap();
    let mut canvas = if term_size.0 > term_size.1 {
        Canvas::new(term_size.1 * 8, term_size.1 * 8, None)
    } else {
        Canvas::new(term_size.0 * 4, term_size.0 * 4, None)
    };

    let mut cam_zoom = settings.cam_zoom;
    let mut cam_xy = 0.;
    let mut cam_z = 0.;

    // set the initial coordinates
    focus_target(settings.coords, 0., &mut cam_xy, &mut cam_z);

    let mut globe = GlobeConfig::new()
        .use_template(GlobeTemplate::Earth)
        .with_camera(CameraConfig::new(cam_zoom, cam_xy, cam_z))
        .display_night(settings.night)
        .build();

    let mut globe_rot_speed = settings.globe_rotation_speed / 1000.;
    let mut cam_rot_speed = settings.cam_rotation_speed / 1000.;

    let mut last_drag_pos = None;
    let mut moving_towards_target: Option<(f32, f32)> = None;

    loop {
        if poll(Duration::from_millis(1000 / settings.refresh_rate as u64)).unwrap() {
            match read().unwrap() {
                Event::Key(event) => match event.code {
                    KeyCode::Char(char) => match char {
                        '-' => globe_rot_speed -= 0.005,
                        '+' => globe_rot_speed += 0.005,
                        ',' => cam_rot_speed -= 0.005,
                        '.' => cam_rot_speed += 0.005,
                        'n' => globe.display_night = !globe.display_night,
                        // vim-style navigation with hjkl
                        'h' => cam_xy += 0.1,
                        'l' => cam_xy -= 0.1,
                        'k' => {
                            if cam_z < 1.5 {
                                cam_z += 0.1;
                            }
                        }
                        'j' => {
                            if cam_z > -1.5 {
                                cam_z -= 0.1;
                            }
                        }
                        _ => break,
                    },
                    KeyCode::PageUp => cam_zoom += 0.1,
                    KeyCode::PageDown => cam_zoom -= 0.1,
                    KeyCode::Up => {
                        if cam_z < 1.5 {
                            cam_z += 0.1;
                        }
                    }
                    KeyCode::Down => {
                        if cam_z > -1.5 {
                            cam_z -= 0.1;
                        }
                    }
                    KeyCode::Left => cam_xy += 0.1,
                    KeyCode::Right => cam_xy -= 0.1,
                    KeyCode::Enter => {
                        focus_target(settings.coords, globe.angle / 2., &mut cam_xy, &mut cam_z);
                        // moving_towards_target = Some(settings.coords);
                    }
                    _ => (),
                },
                Event::Mouse(event) => match event {
                    MouseEvent::Drag(_, x, y, _) => {
                        if let Some(last) = last_drag_pos {
                            let (x_last, y_last) = last;
                            let x_diff = x as Float - x_last as Float;
                            let y_diff = y as Float - y_last as Float;

                            if y_diff > 0. && cam_z < 1.5 {
                                cam_z += 0.1;
                            } else if y_diff < 0. && cam_z > -1.5 {
                                cam_z -= 0.1;
                            }

                            cam_xy += x_diff * PI / 30.;
                            cam_xy += y_diff * PI / 30.;
                        }
                        last_drag_pos = Some((x, y))
                    }
                    MouseEvent::ScrollUp(..) => cam_zoom -= 0.1,
                    MouseEvent::ScrollDown(..) => cam_zoom += 0.1,
                    _ => last_drag_pos = None,
                },
                Event::Resize(width, height) => {
                    term_size = (width, height);
                    canvas = if width > height {
                        Canvas::new(height * 8, height * 8, None)
                    } else {
                        Canvas::new(width * 4, width * 4, None)
                    };
                }
            }
        }

        // apply globe rotation
        globe.angle += globe_rot_speed;
        cam_xy -= globe_rot_speed / 2.;

        // apply camera rotation
        cam_xy -= cam_rot_speed;

        // clip camera zoom
        if cam_zoom < 1.0 {
            cam_zoom = 1.0;
        }

        if let Some(target_coords) = moving_towards_target {
            if move_towards_target(
                settings.focus_speed,
                target_coords,
                cam_zoom,
                globe.angle / 2.,
                &mut cam_xy,
                &mut cam_z,
                &mut cam_zoom,
            ) {
                moving_towards_target = None;
            }
        }

        globe.camera.update(cam_zoom, cam_xy, cam_z);

        // render globe on the canvas
        canvas.clear();
        globe.render_on(&mut canvas);

        // print canvas to terminal
        print_canvas(&mut canvas, &term_size, &mut stdout);
    }

    stdout.execute(cursor::Show).unwrap();
    stdout.execute(cursor::EnableBlinking).unwrap();
    stdout
        .execute(crossterm::event::DisableMouseCapture)
        .unwrap();

    terminal::disable_raw_mode().unwrap();
    stdout.execute(terminal::Clear(ClearType::All)).unwrap();
}

//TODO animate zoom
/// Rotates the camera towards given target coordinates.
pub fn move_towards_target(
    speed: f32,
    coords: (f32, f32),
    target_zoom: f32,
    xy_offset: f32,
    cam_xy: &mut f32,
    cam_z: &mut f32,
    cam_zoom: &mut f32,
) -> bool {
    let (cx, cy) = coords;
    let target_xy = (cx * PI - xy_offset) * -1. - 1.5;
    let target_z = cy * 3. - 1.5;

    let diff_xy = target_xy - *cam_xy;
    let diff_z = target_z - *cam_z;

    if diff_xy.abs() < 0.01 && diff_z.abs() < 0.01 {
        return true;
    }

    let mut xy_move = 0.01 * speed + (diff_xy.abs() / 30. * speed);
    if diff_xy.abs() < 0.07 {
        xy_move = xy_move / 5.;
    }
    if diff_xy > 0. {
        *cam_xy += xy_move;
    } else if diff_xy < 0. {
        *cam_xy -= xy_move;
    }

    let mut z_move = 0.005 * speed + (diff_z.abs() / 30. * speed);
    if diff_z.abs() < 0.07 {
        z_move = z_move / 5.;
    }
    if diff_z > 0. {
        *cam_z += z_move;
    } else if diff_z < 0. {
        *cam_z -= z_move;
    }

    false
}
