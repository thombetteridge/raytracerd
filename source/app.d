import std.stdio;
import std.conv;
import std.math;
import core.memory;

import vec3;

enum float ASPECT_RATIO = 16.0 / 9.0; // Ratio of image width over height
enum int IMAGE_WIDTH = 400; // Rendered image width in pixel count
enum int SAMPLES_PER_PIXEL = 50; // Count of random samples for each pixel
enum int MAX_DEPTH = 20; // Maximum number of ray bounces into scene

enum float VFOV = 30; // Vertical view angle (field of view)
enum Point3 LOOKFROM = Point3(13, 2, 3); // Point camera is looking from
enum Point3 LOOKAT = Point3(0, 0, 0); // Point camera is looking at
enum Vec3 VUP = Vec3(0, 1, 0); // Camera-relative "up" direction

enum float DEFOCUS_ANGLE = 0.6; // Variation angle of rays through each pixel
enum float FOCUS_DIST = 10.0; // Distance from camera LOOKFROM point to plane of perfect focus

enum int IMAGE_HEIGHT = to!int(IMAGE_WIDTH / ASPECT_RATIO); // Rendered image height
enum float PIXEL_SAMPLES_SCALE = 1.0 / SAMPLES_PER_PIXEL; // Color scale factor for a sum of pixel samples
enum Point3 CENTER = LOOKFROM; // Camera center
enum Vec3 PIXEL_DELTA_U = VIEWPORT_U / IMAGE_WIDTH; // Offset to pixel to the right
enum Vec3 PIXEL_DELTA_V = VIEWPORT_V / IMAGE_HEIGHT; // Offset to pixel below
enum Point3 PIXEL00_LOC = VIEWPORT_UPPER_LEFT + (PIXEL_DELTA_U + PIXEL_DELTA_V) * 0.5; // Location of pixel 0, 0
// Camera frame basis vectors
enum Vec3 W = unit_vector(LOOKFROM - LOOKAT);
enum Vec3 U = unit_vector(cross(VUP, W));
enum Vec3 V = cross(W, U); // Calculate the vectors across the horizontal and down the vertical viewport edges.
enum DEFOCUS_RADIUS = FOCUS_DIST * tan(degrees_to_radians(DEFOCUS_ANGLE / 2.0));
enum Vec3 DEFOCUS_DISK_U = U * DEFOCUS_RADIUS; // Defocus disk horizontal radius
enum Vec3 DEFOCUS_DISK_V = V * DEFOCUS_RADIUS; // Defocus disk vertical radius
enum VIEWPORT_UPPER_LEFT = CENTER - (FOCUS_DIST * W) - VIEWPORT_U / 2.0 - VIEWPORT_V / 2.0;

enum VIEWPORT_U = VIEWPORT_WIDTH * U; // Vector across viewport horizontal edge
enum VIEWPORT_V = VIEWPORT_HEIGHT * -V; // Vector down viewport vertical edge

enum THETA = degrees_to_radians(VFOV);
enum H = tan(THETA / 2.0);
enum VIEWPORT_HEIGHT = 2.0 * H * FOCUS_DIST;
enum VIEWPORT_WIDTH = VIEWPORT_HEIGHT * (to!float(IMAGE_WIDTH) / IMAGE_HEIGHT);

void main() {
	GC.disable();

	// World
	Hittable_list world;

	auto ground = Sphere(Point3(0, -1000, 0), 1000);
	ground.mat = MaterialType.lambertian;
	ground.albedo = Colour(0.5, 0.5, 0.5);

	world.add(ground);
	import std.random : uniform01;

	foreach (a; -12 .. 12) {
		foreach (b; -12 .. 12) {
			const choose_mat = uniform01();
			const center = Point3(a + 0.9 * uniform01(), 0.2, b + 0.9 * uniform01());

			if ((center - Point3(4, 0.2, 0)).length > 0.9) {
				if (choose_mat < 0.5) {
					// diffuse
					auto sphere = Sphere(center, 0.2);
					sphere.mat = MaterialType.lambertian;
					sphere.albedo = vec3_random * vec3_random;

					world.add(sphere);
				}
				else if (choose_mat < 0.95) {
					// metal
					auto sphere = Sphere(center, 0.2);
					sphere.mat = MaterialType.metal;
					sphere.albedo = vec3_random(0.5, 1.0);
					//sphere.albedo = Colour(0.7, 0.6, 0.5);

					world.add(sphere);
				}
				else {
					// glass
					auto sphere = Sphere(center, 0.2);
					sphere.mat = MaterialType.dielectric;
					sphere.refraction_index = 1.5;

					world.add(sphere);
				}
			}
		}
	}

	auto sphere_1 = Sphere(Point3(0, 1, 0), 1.0);
	sphere_1.mat = MaterialType.dielectric;
	sphere_1.refraction_index = 1.5;

	auto sphere_2 = Sphere(Point3(-4, 1, 0), 1.0);
	sphere_2.mat = MaterialType.lambertian;
	sphere_2.albedo = Colour(0.4, 0.2, 0.1);

	auto sphere_3 = Sphere(Point3(4, 1, 0), 1.0);
	sphere_3.mat = MaterialType.metal;
	sphere_3.albedo = Colour(0.7, 0.6, 0.5);

	world.add(sphere_1);
	world.add(sphere_2);
	world.add(sphere_3);

	Camera cam;
	cam.render(world);
}

float linear_to_gamma(float linear_component) {
	if (linear_component < 0) {
		return 0;
	}
	return linear_component.sqrt;
}

void write_bmp(string filename, Colour[][] pixels) {
	import std.bitmanip;
	import std.file;
	import std.algorithm : clamp, min;
	import std.conv;

	immutable width = pixels[0].length;
	immutable height = pixels.length;

	immutable row_padding = (4 - (width * 3) % 4) % 4;
	immutable row_size = width * 3 + row_padding;
	immutable pixel_data_size = row_size * height;
	immutable file_size = 54 + pixel_data_size;

	// BMP Header (14 bytes)
	ubyte[14] header;
	//header.length = 14;
	header[0 .. 2] = cast(ubyte[])("BM"); // Signature
	header[2 .. 6] = nativeToLittleEndian(to!uint(file_size)); // File size
	header[10 .. 14] = nativeToLittleEndian(to!uint(54)); // Pixel data offset

	// DIB Header (40 bytes)
	ubyte[40] dib;
	//dib.length = 40;
	dib[0 .. 4] = nativeToLittleEndian(to!uint(40)); // DIB header size
	dib[4 .. 8] = nativeToLittleEndian(to!int(width));
	dib[8 .. 12] = nativeToLittleEndian(to!int(height));
	dib[12 .. 14] = nativeToLittleEndian(to!ushort(1)); // Color planes
	dib[14 .. 16] = nativeToLittleEndian(to!ushort(24)); // Bits per pixel
	dib[20 .. 24] = nativeToLittleEndian(to!uint(pixel_data_size)); // Image size

	// Build pixel data (BMP stores bottom to top)
	scope pixel_data = new ubyte[](pixel_data_size);

	foreach (y; 0 .. height) {
		const row = pixels[height - 1 - y]; // BMP is bottom-up
		foreach (x; 0 .. width) {
			enum intensity = Interval(0.000, 0.999);

			const c = row[x];
			const r = linear_to_gamma(c.r);
			const g = linear_to_gamma(c.g);
			const b = linear_to_gamma(c.b);

			const rbyte = to!ubyte(min(255,255.999 * intensity.clamp(r)));
			const gbyte = to!ubyte(min(255,255.999 * intensity.clamp(g)));
			const bbyte = to!ubyte(min(255,255.999 * intensity.clamp(b)));

			immutable i = y * row_size + x * 3;
			pixel_data[i + 0] = bbyte;
			pixel_data[i + 1] = gbyte;
			pixel_data[i + 2] = rbyte;
		}

		// padding
		foreach (p; 0 .. row_padding) {
			pixel_data[y * row_size + width * 3 + p] = 0;
		}
	}

	// Write to file
	const bmp = header ~ dib ~ pixel_data;
	write(filename, bmp);
}

void write_colour(in Colour pixel_colour) {
	static const intensity = Interval(0.000, 0.999);
	const r = linear_to_gamma(pixel_colour.r);
	const g = linear_to_gamma(pixel_colour.g);
	const b = linear_to_gamma(pixel_colour.b);

	const rbyte = to!int(255.999 * intensity.clamp(r));
	const gbyte = to!int(255.999 * intensity.clamp(g));
	const bbyte = to!int(255.999 * intensity.clamp(b));
	writeln(rbyte, ' ', gbyte, ' ', bbyte);
}

// Ray
struct Ray {
	Point3 origin;
	Vec3 direction;
	this(in Point3 origin, in Vec3 direction) {
		this.origin = Point3(origin);
		this.direction = Vec3(direction);
	}

	Point3 at(float t) const {
		return origin + (direction * t);
	}
}

// Hit_record
struct Hit_record {
	Point3 p;
	Vec3 normal;
	Colour albedo;
	float refraction_index;
	float t;
	MaterialType mat;
	bool front_face;
	bool valid;

	void set_face_normal(in Ray r, in Vec3 outward_normal) {
		front_face = dot(r.direction, outward_normal) < 0;
		normal = front_face ? outward_normal : -outward_normal;
		valid = true;
	}
}

// Shape Structs

struct Sphere {
	Point3 center;
	Colour albedo;
	float radius;
	float refraction_index;
	MaterialType mat;

	this(in Point3 center, float radius) {
		import std.algorithm : max;

		this.center = center;
		this.radius = max(0, radius);
	}

	Hit_record hit(in Ray r, in Interval ray_t) const {
		const oc = center - r.origin;
		const a = r.direction.length_squared;
		const H = dot(r.direction, oc);
		const c = oc.length_squared - radius * radius;

		const discriminant = H * H - a * c;
		if (discriminant < 0) {
			return Hit_record(); // valid = false
		}

		const sqrtd = sqrt(discriminant);
		float root = (H - sqrtd) / a;
		if (!ray_t.surrounds(root)) {
			root = (H + sqrtd) / a;
			if (!ray_t.surrounds(root))
				return Hit_record();
		}

		Hit_record rec;
		rec.t = root;
		rec.p = r.at(rec.t);
		const outward_normal = (rec.p - center) / radius;
		rec.set_face_normal(r, outward_normal);
		rec.mat = mat;
		rec.albedo = albedo;
		rec.refraction_index = refraction_index;
		return rec;
	}
}

// Hittable Template
struct Hittable(T) {
	T shape;
	this(T s) {
		shape = s;
	}

	Hit_record hit(in Ray r, in Interval ray_t) const {
		return shape.hit(r, ray_t);
	}
}

// Hittable_list
struct Hittable_list {
	Hittable!Sphere[] spheres;

	void add(T)(T object) {
		static if (is(T == Sphere)) {
			spheres ~= Hittable!Sphere(object);
		}
		else {
			static assert(false, "Unsupported shape type: " ~ T.stringof);
		}
	}

	void clear() {
		spheres.length = 0;
	}

	Hit_record hit(in Ray r, in Interval ray_t) const {
		Hit_record closest_rec;
		Interval current_t = ray_t;

		foreach (ref s; spheres) {
			auto temp_rec = s.hit(r, current_t);
			if (temp_rec.valid) {
				current_t.max = temp_rec.t;
				closest_rec = temp_rec;
			}
		}
		return closest_rec;
	}
}

enum MaterialType : ubyte {
	lambertian,
	metal,
	dielectric
}

bool scatter(in Ray r_in, in Hit_record rec, out Colour attenuation, out Ray scattered) {
	switch (rec.mat) with (MaterialType) {
	case lambertian:
		return lambertian_scatter(r_in, rec, attenuation, scattered);
	case metal:
		return metal_scatter(r_in, rec, attenuation, scattered);
	case dielectric:
		return dielectric_scatter(r_in, rec, attenuation, scattered);
	default:
		return false;
	}
}

bool lambertian_scatter(in Ray r_in, in Hit_record rec, out Colour attenuation, out Ray scattered) {
	auto scatter_direction = rec.normal + random_unit_vector();
	if (scatter_direction.near_zero) {
		scatter_direction = rec.normal;
	}
	scattered = Ray(rec.p, scatter_direction);
	attenuation = rec.albedo;
	return true;
}

bool metal_scatter(in Ray r_in, in Hit_record rec, out Colour attenuation, out Ray scattered) {
	const reflected = reflect(r_in.direction, rec.normal);
	scattered = Ray(rec.p, reflected);
	attenuation = rec.albedo;
	return true;
}

bool dielectric_scatter(in Ray r_in, in Hit_record rec, out Colour attenuation, out Ray scattered) {
	import std.algorithm : min;
	import std.random : uniform01;

	const refraction_index = rec.refraction_index;

	attenuation = Colour(1.0, 1.0, 1.0);

	const ri = rec.front_face ? (1.0 / refraction_index) : refraction_index;
	const unit_direction = unit_vector(r_in.direction);
	const cos_theta = min(dot(-unit_direction, rec.normal), 1.0);
	const sin_theta = sqrt(1.0 - cos_theta * cos_theta);
	const cannot_refract = ri * sin_theta > 1.0;

	Vec3 direction;
	if (cannot_refract || reflectance(cos_theta, ri) > uniform01!float) {
		direction = reflect(unit_direction, rec.normal);
	}
	else {
		direction = refract(unit_direction, rec.normal, ri);
	}
	scattered = Ray(rec.p, direction);
	return true;
}

float reflectance(float cosine, float refraction_index) {
	import std.math : pow;

	// Use Schlick's approximation for reflectance.
	auto r0 = (1 - refraction_index) / (1 + refraction_index);

	r0 = r0 * r0;
	return r0 + (1 - r0) * pow((1 - cosine), 5);
}

struct Interval {
	float min = float.infinity;
	float max = -float.infinity;

	this(float min, float max) {
		this.min = min;
		this.max = max;
	}

	auto size() const => max - min;
	auto contains(float x) const => min <= x && x <= max;
	auto surrounds(float x) const => min < x && x < max;
	
	float clamp(float x) const {
		if (x < min)
			return min;
		if (x > max)
			return max;
		return x;
	}
}

float degrees_to_radians(float degrees) {
	enum pi_180 = to!float(PI / 180.0f);
	return degrees * pi_180;
}

const empty = Interval(+float.infinity, -float.infinity);
const universe = Interval(-float.infinity, +float.infinity);

struct Camera {

	void render(in Hittable_list world) {
		import std.range : iota;
		import std.algorithm : fold;
		import std.parallelism : parallel;
		import core.thread;
		import core.atomic;

		shared int doneRows;
		immutable totalRows = IMAGE_HEIGHT;

		auto progressThread = new Thread({
			while (true) {
				int current = atomicLoad(doneRows);
				stderr.writef("\rProgress: %d%%", std.conv.to!int(
					(std.conv.to!float(current) / totalRows) * 100));
				stderr.flush();
				if (current >= totalRows) {
					break;
				}
				Thread.sleep(500.msecs);
			}
		});
		progressThread.start();

		auto pixels = new Colour[][](IMAGE_HEIGHT);
		foreach (j; iota(IMAGE_HEIGHT).parallel) {
			Colour[] row = new Colour[](IMAGE_WIDTH);
			foreach (i; iota(IMAGE_WIDTH).parallel) {
				const pixel_colour = iota(SAMPLES_PER_PIXEL)
					.fold!((result, e) => result + ray_color(get_ray(i, j), MAX_DEPTH, world))(
						Colour(0, 0, 0));

				row[i] = pixel_colour * PIXEL_SAMPLES_SCALE;
			}
			pixels[j] = row;
			atomicOp!"+="(doneRows, 1);
		}

		//write_out_pixels(pixels);
		write_bmp("output.bmp", pixels); // Save BMP
		progressThread.join();
		stderr.write("\rDone.                   \n");
		stderr.flush();
	}

	void write_out_pixels(in Colour[][] pixels) {
		writef("P3\n%s %s\n255\n", IMAGE_WIDTH, IMAGE_HEIGHT);
		foreach (row; pixels) {
			foreach (pixel; row) {
				write_colour(pixel);
			}
		}
	}

	Ray get_ray(int i, int j) const {
		// Construct a camera ray originating from the origin and directed at randomly sampled
		// point around the pixel location i, j.

		const offset = sample_square();
		const pixel_sample = PIXEL00_LOC + ((i + offset.x) * PIXEL_DELTA_U) + (
			(j + offset.y) * PIXEL_DELTA_V);

		const ray_origin = (DEFOCUS_ANGLE <= 0) ? CENTER : defocus_disk_sample();
		const ray_direction = pixel_sample - ray_origin;

		return Ray(ray_origin, ray_direction);
	}

	Point3 defocus_disk_sample() const {
		// Returns a random point in the camera defocus disk.
		auto p = random_in_unit_disk();
		return CENTER + (p.x * DEFOCUS_DISK_U) + (p.y * DEFOCUS_DISK_V);
	}

	Vec3 sample_square() const {
		import std.random : uniform01;

		// Returns the vector to a random point in the [-.5,-.5]-[+.5,+.5] unit square.
		return Vec3(uniform01 - 0.5, uniform01 - 0.5, 0);
	}

	Colour ray_color(in Ray r, int depth, in Hittable_list world) {
		if (depth == 0) {
			return Colour(0, 0, 0);
		}
		auto rec = world.hit(r, Interval(0.001, float.infinity));
		if (rec.valid) {
			Ray scattered;
			Colour attenuation;

			if (scatter(r, rec, attenuation, scattered)) {
				return attenuation * ray_color(scattered, depth - 1, world);
			}
			return Colour(0, 0, 0);
		}
		const unit_direction = unit_vector(r.direction);
		const a = 0.5 * (unit_direction.y + 1.0);
		return Colour(1.0, 1.0, 1.0) * (1.0 - a) + Colour(0.5, 0.7, 1.0) * a;
	}
}
