// Vec3, Point3, Colour
struct Vec3
{
	double x, y, z;

	this(double x_, double y_, double z_)
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

	double r() const
	{
		return x;
	}

	double g() const
	{
		return y;
	}

	double b() const
	{
		return z;
	}

	string toString() const
	{
		import std.format : format;

		return format("%s %s %s", x, y, z);
	}

	ref Vec3 opOpAssign(string op)(Vec3 rhs)
			if (op == "+" || op == "-" || op == "*" || op == "/")
	{
		mixin("x ", op, "= rhs.x");
		mixin("y ", op, "= rhs.y");
		mixin("z ", op, "= rhs.z");
		return this;
	}

	ref Vec3 opOpAssign(string op)(double scalar)
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

	Vec3 opBinary(string op)(double scalar) const
	if (op == "+" || op == "-" || op == "*" || op == "/")
	{
		return mixin("Vec3(scalar ", op, " x, scalar ", op, " y, scalar ", op, " z)");
	}

	Vec3 opBinaryRight(string op : "*")(double scalar) const
	if (op == "+" || op == "-" || op == "*" || op == "/")
	{
		return mixin("Vec3(x ", op, " scalar, y ", op, "  scalar, z ", op, " scalar)");
	}

	double length() const
	{
		import std.math : sqrt;

		return sqrt(length_squared);
	}

	double length_squared() const
	{
		return x * x + y * y + z * z;
	}

	bool near_zero() const
	{
		import std.math : abs;

		// Return true if the vector is close to zero in all dimensions.
		static const double s = 1e-8;
		return (x.abs < s) && (y.abs < s) && (z.abs < s);
	}
}

Vec3 vec3_random()
{
	import std.random : uniform01;

	return Vec3(uniform01, uniform01, uniform01);
}

Vec3 vec3_random(double min, double max)
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

double dot(Vec3 a, Vec3 b)
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

Vec3 refract(in Vec3 a, in Vec3 b, double etai_over_etat)
{
	import std.algorithm : min;
	import std.math : sqrt, abs;

	const cos_theta = min(dot(-a, b), 1.0);
	const r_out_perp = etai_over_etat * (
		a + cos_theta * b);
	const r_out_parallel = -sqrt(
		(1.0 - r_out_perp.length_squared).abs) * b;
	return r_out_perp + r_out_parallel;
}

Vec3 random_in_unit_disk()
{
	while (true)
	{
		import std.random : uniform;

		const p = Vec3(uniform(-1.0, 1.0), uniform(-1.0, 1.0), 0);
		if (p.length_squared() < 1)
			return p;
	}
}
