import std.stdio;
import std.conv;
import std.math;

import core.memory;

import vec3;

void main()
{
	GC.disable();

	// World
	Hittable_list world;

	auto ground = Sphere(Point3(0, -1000, 0), 1000);
	ground.mat = MaterialType.lambertian;
	ground.albedo = Colour(0.5, 0.5, 0.5);

	world.add(ground);

	import std.random : uniform01;

	foreach (a; -12 .. 12)
	{
		foreach (b; -12 .. 12)
		{
			const choose_mat = uniform01();
			const center = Point3(a + 0.9 * uniform01(), 0.2, b + 0.9 * uniform01());

			if ((center - Point3(4, 0.2, 0)).length > 0.9)
			{
				if (choose_mat < 0.8)
				{
					// diffuse
					auto sphere = Sphere(center, 0.2);
					sphere.mat = MaterialType.lambertian;
					sphere.albedo = vec3_random * vec3_random;

					world.add(sphere);
				}
				else if (choose_mat < 0.95)
				{
					// metal
					auto sphere = Sphere(center, 0.2);
					sphere.mat = MaterialType.metal;
					sphere.albedo = vec3_random(0.5, 1);

					world.add(sphere);
				}
				else
				{
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
	cam.aspect_ratio = 16.0 / 9.0;
	cam.image_width = 1600;
	cam.samples_per_pixel = 100;
	cam.max_depth = 50;
	cam.vfov = 20;
	cam.lookfrom = Point3(13, 2, 3);
	cam.lookat = Point3(0, 0, 0);
	cam.vup = Vec3(0, 1, 0);
	cam.defocus_angle = 0.6;
	cam.focus_dist = 10.0;
	cam.render(world);
}

float linear_to_gamma(float linear_component)
{
	if (linear_component < 0)
	{
		return 0;
	}
	return linear_component.sqrt;
}

void write_colour(in Colour pixel_colour)
{
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
struct Ray
{
	Point3 origin;
	Vec3 direction;
	this(in Point3 origin, in Vec3 direction)
	{
		this.origin = Point3(origin);
		this.direction = Vec3(direction);
	}

	Point3 at(float t) const
	{
		return origin + (direction * t);
	}
}

// Hit_record
struct Hit_record
{
	Point3 p;
	Vec3 normal;
	Colour albedo;
	float refraction_index;
	float t;
	MaterialType mat;
	bool front_face;
	bool valid;

	void set_face_normal(in Ray r, in Vec3 outward_normal)
	{
		front_face = dot(r.direction, outward_normal) < 0;
		normal = front_face ? outward_normal : -outward_normal;
		valid = true;
	}
}

// Shape Structs

struct Sphere
{
	Point3 center;
	Colour albedo;
	float radius;
	float refraction_index;
	MaterialType mat;

	this(in Point3 center, float radius)
	{
		import std.algorithm : max;

		this.center = center;
		this.radius = max(0, radius);
	}

	Hit_record hit(in Ray r, in Interval ray_t) const
	{
		const oc = center - r.origin;
		const a = r.direction.length_squared;
		const h = dot(r.direction, oc);
		const c = oc.length_squared - radius * radius;

		const discriminant = h * h - a * c;
		if (discriminant < 0)
		{
			return Hit_record(); // valid = false
		}

		const sqrtd = sqrt(discriminant);
		float root = (h - sqrtd) / a;
		if (!ray_t.surrounds(root))
		{
			root = (h + sqrtd) / a;
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
struct Hittable(T)
{
	T shape;
	this(T s)
	{
		shape = s;
	}

	Hit_record hit(in Ray r, in Interval ray_t) const
	{
		return shape.hit(r, ray_t);
	}
}

// Hittable_list
struct Hittable_list
{
	Hittable!Sphere[] spheres;

	void add(T)(T object)
	{
		static if (is(T == Sphere))
		{
			spheres ~= Hittable!Sphere(object);
		}
		else
		{
			static assert(false, "Unsupported shape type: " ~ T.stringof);
		}
	}

	void clear()
	{
		spheres.length = 0;
	}

	Hit_record hit(in Ray r, in Interval ray_t) const
	{
		Hit_record closest_rec;
		Interval current_t = ray_t;

		foreach (ref s; spheres)
		{
			auto temp_rec = s.hit(r, current_t);
			if (temp_rec.valid)
			{
				current_t.max = temp_rec.t;
				closest_rec = temp_rec;
			}
		}
		return closest_rec;
	}
}

enum MaterialType : ubyte
{
	lambertian,
	metal,
	dielectric
}

bool scatter(in Ray r_in, in Hit_record rec, out Colour attenuation, out Ray scattered)
{
	if (rec.mat == MaterialType.lambertian)
	{
		return lambertian_scatter(r_in, rec, attenuation, scattered);
	}
	else if (rec.mat == MaterialType.metal)
	{
		return metal_scatter(r_in, rec, attenuation, scattered);
	}
	else if (rec.mat == MaterialType.dielectric)
	{
		return dielectric_scatter(r_in, rec, attenuation, scattered);
	}
	return false;
}

bool lambertian_scatter(in Ray r_in, in Hit_record rec, out Colour attenuation, out Ray scattered)
{
	auto scatter_direction = rec.normal + random_unit_vector();
	if (scatter_direction.near_zero)
	{
		scatter_direction = rec.normal;
	}
	scattered = Ray(rec.p, scatter_direction);
	attenuation = rec.albedo;
	return true;
}

bool metal_scatter(in Ray r_in, in Hit_record rec, out Colour attenuation, out Ray scattered)
{
	const reflected = reflect(r_in.direction, rec.normal);
	scattered = Ray(rec.p, reflected);
	attenuation = rec.albedo;
	return true;
}

bool dielectric_scatter(in Ray r_in, in Hit_record rec, out Colour attenuation, out Ray scattered)
{
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
	if (cannot_refract || reflectance(cos_theta, ri) > uniform01!float)
	{
		direction = reflect(unit_direction, rec.normal);
	}
	else
	{
		direction = refract(unit_direction, rec.normal, ri);
	}
	scattered = Ray(rec.p, direction);
	return true;
}

float reflectance(float cosine, float refraction_index)
{
	import std.math : pow;

	// Use Schlick's approximation for reflectance.
	auto r0 = (1 - refraction_index) / (1 + refraction_index);

	r0 = r0 * r0;
	return r0 + (1 - r0) * pow((1 - cosine), 5);
}

struct Interval
{
	float min = float.infinity;
	float max = -float.infinity;

	this(float min, float max)
	{
		this.min = min;
		this.max = max;
	}

	float size() const
	{
		return max - min;
	}

	bool contains(float x) const
	{
		return min <= x && x <= max;
	}

	bool surrounds(float x) const
	{
		return min < x && x < max;
	}

	float clamp(float x) const
	{
		if (x < min)
			return min;
		if (x > max)
			return max;
		return x;
	}
}

float degrees_to_radians(float degrees)
{
	//import std.math : pi;

	static pi_180 = PI / 180.0;
	return degrees * pi_180;
}

const empty = Interval(+float.infinity, -float.infinity);
const universe = Interval(-float.infinity, +float.infinity);

struct Camera
{
	float aspect_ratio = 1.0; // Ratio of image width over height
	int image_width = 100; // Rendered image width in pixel count
	int samples_per_pixel = 10; // Count of random samples for each pixel
	int max_depth = 10; // Maximum number of ray bounces into scene

	float vfov = 90; // Vertical view angle (field of view)
	Point3 lookfrom = Point3(0, 0, 0); // Point camera is looking from
	Point3 lookat = Point3(0, 0, -1); // Point camera is looking at
	Vec3 vup = Vec3(0, 1, 0); // Camera-relative "up" direction

	float defocus_angle = 0; // Variation angle of rays through each pixel
	float focus_dist = 10; // Distance from camera lookfrom point to plane of perfect focus

	int image_height; // Rendered image height
	float pixel_samples_scale; // Color scale factor for a sum of pixel samples
	Point3 center; // Camera center
	Point3 pixel00_loc; // Location of pixel 0, 0
	Vec3 pixel_delta_u; // Offset to pixel to the right
	Vec3 pixel_delta_v; // Offset to pixel below
	Vec3 u, v, w; // Camera frame basis vectors
	Vec3 defocus_disk_u; // Defocus disk horizontal radius
	Vec3 defocus_disk_v; // Defocus disk vertical radius

	void initialize()
	{
		image_height = to!int(image_width / aspect_ratio);
		pixel_samples_scale = 1.0 / samples_per_pixel;
		center = lookfrom;

		const theta = degrees_to_radians(vfov);
		const h = tan(theta / 2.0);
		const viewport_height = 2.0 * h * focus_dist;
		const viewport_width = viewport_height * (to!float(image_width) / image_height);

		// Calculate the u,v,w unit basis vectors for the camera coordinate frame.
		w = unit_vector(lookfrom - lookat);
		u = unit_vector(cross(vup, w));
		v = cross(w, u); // Calculate the vectors across the horizontal and down the vertical viewport edges.
		const viewport_u = viewport_width * u; // Vector across viewport horizontal edge
		const viewport_v = viewport_height * -v; // Vector down viewport vertical edge

		pixel_delta_u = viewport_u / image_width;
		pixel_delta_v = viewport_v / image_height;

		// Calculate the location of the upper left pixel.
		const viewport_upper_left = center - (focus_dist * w) - viewport_u / 2.0 - viewport_v / 2.0;
		pixel00_loc = viewport_upper_left + (pixel_delta_u + pixel_delta_v) * 0.5;

		// Calculate the camera defocus disk basis vectors.
		const defocus_radius = focus_dist * tan(degrees_to_radians(defocus_angle / 2.0));
		defocus_disk_u = u * defocus_radius;
		defocus_disk_v = v * defocus_radius;
	}

	void render(in Hittable_list world)
	{
		import std.range : iota;
		import std.algorithm : fold;
		import std.parallelism : parallel;
		import core.thread;
		import core.atomic;

		initialize();

		shared int doneRows;
		immutable totalRows = image_height;

		auto progressThread = new Thread({
			while (true)
			{
				int current = atomicLoad(doneRows);
				stderr.writef("\rProgress: %d%%", cast(int)((cast(float) current / totalRows) * 100));
				stderr.flush();
				if (current >= totalRows)
				{
					break;
				}
				Thread.sleep(500.msecs);
			}
		});
		progressThread.start();

		auto pixels = new Colour[][](image_height);
		foreach (j; iota(image_height).parallel)
		{
			Colour[] row = new Colour[](image_width);
			const local_max_depth = max_depth;
			const local_world = world;
			foreach (i; 0 .. image_width)
			{
				const pixel_colour = iota(samples_per_pixel)
					.fold!((result, e) => result + ray_color(get_ray(i, j), local_max_depth, local_world))(
						Colour(0, 0, 0));

				row[i] = pixel_colour * pixel_samples_scale;
			}
			pixels[j] = row;
			atomicOp!"+="(doneRows, 1);
		}

		write_out_pixels(pixels);

		progressThread.join();
		stderr.write("\rDone.                   \n");
		stderr.flush();
	}

	void write_out_pixels(in Colour[][] pixels)
	{
		writef("P3\n%s %s\n255\n", image_width, image_height);
		foreach (row; pixels)
		{
			foreach (pixel; row)
			{
				write_colour(pixel);
			}
		}
	}

	Ray get_ray(int i, int j) const
	{
		// Construct a camera ray originating from the origin and directed at randomly sampled
		// point around the pixel location i, j.

		const offset = sample_square();
		const pixel_sample = pixel00_loc + ((i + offset.x) * pixel_delta_u) + (
			(j + offset.y) * pixel_delta_v);

		const ray_origin = (defocus_angle <= 0) ? center : defocus_disk_sample();
		const ray_direction = pixel_sample - ray_origin;

		return Ray(ray_origin, ray_direction);
	}

	Point3 defocus_disk_sample() const
	{
		// Returns a random point in the camera defocus disk.
		auto p = random_in_unit_disk();
		return center + (p.x * defocus_disk_u) + (p.y * defocus_disk_v);
	}

	Vec3 sample_square() const
	{
		import std.random : uniform01;

		// Returns the vector to a random point in the [-.5,-.5]-[+.5,+.5] unit square.
		return Vec3(uniform01 - 0.5, uniform01 - 0.5, 0);
	}

	Colour ray_color(in Ray r, int depth, in Hittable_list world)
	{
		if (depth == 0)
		{
			return Colour(0, 0, 0);
		}
		auto rec = world.hit(r, Interval(0.001, float.infinity));
		if (rec.valid)
		{
			Ray scattered;
			Colour attenuation;

			if (scatter(r, rec, attenuation, scattered))
			{
				return attenuation * ray_color(scattered, depth - 1, world);
			}
			return Colour(0, 0, 0);
		}
		const unit_direction = unit_vector(r.direction);
		const a = 0.5 * (unit_direction.y + 1.0);
		return Colour(1.0, 1.0, 1.0) * (1.0 - a) + Colour(0.5, 0.7, 1.0) * a;
	}
}
