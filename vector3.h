#ifndef H_VECTOR3
#define H_VECTOR3

class Vector3 {
 private:
  double coords[3];

 public:
  double x;
  double y;
  double z;
  Vector3();
  Vector3(double x, double y, double z);
  Vector3& operator=(const Vector3& rhs);
  ~Vector3();

  void set_coords(double x, double y, double z);
  void return_coords(double & x, double & y, double & z);
  // void make_unitvector();
  double dot(const Vector3 & vec);
};



#endif
