// Vec3, Point3, Colour
struct Vec3
{
	float x, y, z;

	this(float x_, float y_, float z_)
	{
		this.x = x_;
		this.y = y_;
		this.z = z_;
	}

	this(in Vec3 rhs)
	{
		this.x = rhs.x;
		this.y = rhs.y;
		this.z = rhs.z;
	}

	auto ref r()
	{
		return x;
	}

	auto ref g()
	{
		return y;
	}

	auto ref b()
	{
		return z;
	}

	auto r() const
	{
		return x;
	}

	auto g() const
	{
		return y;
	}

	auto b() const
	{
		return z;
	}

	string toString() const
	{
		import std.format : format;

		return format("%s %s %s", x, y, z);
	}

	ref Vec3 opOpAssign(string op)(in Vec3 rhs)
			if (op == "+" || op == "-" || op == "*" || op == "/")
	{
		mixin("x ", op, "= rhs.x");
		mixin("y ", op, "= rhs.y");
		mixin("z ", op, "= rhs.z");
		return this;
	}

	ref Vec3 opOpAssign(string op)(float scalar)
			if (op == "+" || op == "-" || op == "*" || op == "/")
	{
		mixin("x ", op, "= scalar");
		mixin("y ", op, "= scalar");
		mixin("z ", op, "= scalar");
		return this;
	}

	Vec3 opUnary(string op : "-")() const
	{
		return Vec3(-x, -y, -z);
	}

	Vec3 opBinary(string op)(Vec3 rhs) const
	if (op == "+" || op == "-" || op == "*" || op == "/")
	{
		return mixin(" Vec3( this.x ", op, " rhs.x, this.y ", op, " rhs.y, this.z ", op, " rhs.z)");
	}

	Vec3 opBinary(string op)(float scalar) const
	if (op == "+" || op == "-" || op == "*" || op == "/")
	{
		return mixin("Vec3(x ", op, " scalar, y ", op, " scalar, z ", op, " scalar)");
	}

	Vec3 opBinaryRight(string op : "*")(float scalar) const
	if (op == "+" || op == "-" || op == "*" || op == "/")
	{
		return this.opBinary!op(scalar);
	}

	auto length() const
	{
		import std.math : sqrt;

		return sqrt(length_squared);
	}

	auto length_squared() const
	{
		return x * x + y * y + z * z;
	}

	bool near_zero() const
	{
		import std.math : abs;

		// Return true if the vector is close to zero in all dimensions.
		static const float s = 1e-8;
		return (x.abs < s) && (y.abs < s) && (z.abs < s);
	}
}

Vec3 vec3_random()
{
	import std.random : uniform01;

	return Vec3(uniform01, uniform01, uniform01);
}

Vec3 vec3_random(float min, float max)
{
	import std.random : uniform;

	return Vec3(uniform(min, max), uniform(min, max), uniform(min, max));
}

alias Point3 = Vec3;
alias Colour = Vec3;

Vec3 unit_vector(Vec3 v)
{
	return v / v.length;
}

Vec3 random_unit_vector()
{
	import std.math : sqrt;

	while (true)
	{
		const p = vec3_random(-1, 1);
		const lensq = p.length_squared();
		if (1e-160 < lensq && lensq <= 1)
			return p / sqrt(lensq);
	}
}

Vec3 random_on_hemisphere(in Vec3 normal)
{
	auto on_unit_sphere = random_unit_vector();
	if (dot(on_unit_sphere, normal) > 0.0) // In the same hemisphere as the normal
		return on_unit_sphere;
	else
		return -on_unit_sphere;
}

auto dot(in Vec3 a, in Vec3 b) pure
{
	return a.x * b.x + a.y * b.y + a.z * b.z;
}

Vec3 cross(in Vec3 a, in Vec3 b)
{
	return Vec3(
		a.y * b.z - a.z * b.y,
		a.z * b.x - a.x * b.z,
		a.x * b.y - a.y * b.x);
}

Vec3 reflect(in Vec3 a, in Vec3 b)
{
	return a - 2 * dot(a, b) * b;
}

Vec3 refract(in Vec3 a, in Vec3 b, float etai_over_etat)
{
	import std.algorithm : min;
	import std.math : sqrt, abs;

	const cos_theta = min(dot(-a, b), 1.0);
	const r_out_perp = etai_over_etat * (a + cos_theta * b);
	const r_out_parallel = -sqrt((1.0 - r_out_perp.length_squared).abs) * b;
	return r_out_perp + r_out_parallel;
}

Vec3 random_in_unit_disk()
{
	while (true)
	{
		import std.random : uniform;

		const p = Vec3(uniform(-1.0, 1.0), uniform(-1.0, 1.0), 0);
		if (p.length_squared() < 1)
		{
			return p;
		}
	}
}

unittest
{
	import std.math : abs, approxEqual;

	Vec3 a = Vec3(1, 2, 3);
	Vec3 b = Vec3(4, 5, 6);

	// Constructor & accessors
	assert(a.r == 1);
	assert(a.g == 2);
	assert(a.b == 3);

	// Copy constructor
	const c = Vec3(a);
	assert(c.x == a.x && c.y == a.y && c.z == a.z);

	// Operator overloads (vec3)
	assert((a + b).toString() == "5 7 9");
	assert((b - a).toString() == "3 3 3");
	assert((a * b).toString() == "4 10 18");
	assert((b / a).toString() == "4 2.5 2");

	// Operator overloads (scalar)
	assert((a * 2).toString() == "2 4 6");
	assert((2 * a).toString() == "2 4 6");
	assert((a / 2).toString() == "0.5 1 1.5");

	// Unary
	assert((-a).toString() == "-1 -2 -3");

	// Length
	assert(approxEqual(a.length, 3.74165, 1e-5));
	assert(a.length_squared == 14);

	// near_zero
	assert(!a.near_zero);
	assert(Vec3(1e-9, 0, 0).near_zero);

	// dot
	assert(dot(a, b) == 32);

	// cross
	const cr = cross(a, b);
	assert(cr.toString() == "-3 6 -3");

	// reflect
	const r = reflect(Vec3(1, -1, 0), Vec3(0, 1, 0));
	assert(r.toString() == "1 1 0");

	// refract sanity check
	const refracted = refract(Vec3(1, 1, 0), Vec3(0, 1, 0), 1.0 / 1.5);
	assert(refracted.length <= 1); // should not exceed unit length

	// unit_vector
	const u = unit_vector(Vec3(3, 0, 0));
	assert(u.toString() == "1 0 0");

	// random_unit_vector
	foreach (i; 0 .. 5)
	{
		const v = random_unit_vector();
		assert(approxEqual(v.length, 1, 1e-8));
	}

	// random_on_hemisphere
	foreach (i; 0 .. 5)
	{
		const n = Vec3(0, 1, 0);
		const v = random_on_hemisphere(n);
		assert(dot(v, n) >= 0);
	}

	// random_in_unit_disk
	foreach (i; 0 .. 5)
	{
		const v = random_in_unit_disk();
		assert(v.z == 0);
		assert(v.length_squared < 1);
	}
}
