// tracer.cpp : Defines the entry point for the console application.
//
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <windows.h>
#include <algorithm>

using namespace std;

#include "vector3d.h"
#include "colors.h"
#include "MATRIX3D.H"
#include "matrix.h"

const bool AA=true;

const Vector3D colorIntensity=Vector3D(0.299,0.578,0.114);


struct Ray {
	Vector3D origin;
	Vector3D dir;
};

class Camera {
public:
	Vector3D location;
	Vector3D forward,right,up;
	float width,height;
	Camera(Vector3D loc,Vector3D lookAt,float w,float h,Vector3D u=Vector3D(0,0,1)) {
		location=loc;
		forward=lookAt-location; forward.normalize();
		right=forward^u; right.normalize();
		up=right^forward; up.normalize();
		width=w; 
		height=h;
	}
	Ray getRay(float x,float y) {
		Ray r;
		r.origin=location;
		r.dir=forward+(x-0.5)*width*right+(0.5-y)*height*up;
		r.dir.normalize();
		return r;
	}
};

Camera *camera;

struct Surface {
	Vector3D n;
	Vector3D ambient;
	Vector3D diffuse;
	Vector3D specular;
	float phong;
	Vector3D reflective;
	Vector3D refractive;
	float nRefr;
};

class Material {
public:
	virtual Surface getSurface(Vector3D pos,Vector3D n)=0;
};

class UniformMaterial : public Material {
private:
	Vector3D ambient, diffuse, specular, reflective, refractive;
	float phong;
	float nRefr;
public:
	UniformMaterial(Vector3D a,Vector3D d,Vector3D s=Black,float p=0,Vector3D r=Black,Vector3D refr=Black,float nr=1) {
		ambient=a; diffuse=d; specular=s; phong=p; reflective=r;
		refractive=refr; nRefr=nr;
	}
	virtual Surface getSurface(Vector3D pos,Vector3D n) {
		Surface s;
		s.n=n;
		s.ambient=ambient;
		s.diffuse=diffuse;
		s.specular=specular;
		s.phong=phong;
		s.reflective=reflective;
		s.refractive=refractive;
		s.nRefr=nRefr;
		return s;
	}
};

class Checkers : public Material {
private:
	Material *black,*white;
	float size;
public:
	Checkers(Material *b,Material *w,float s) {
		black=b; white=w; size=s;
	}
	virtual Surface getSurface(Vector3D pos,Vector3D n) {
		pos/=size;
		if (pos.x<0) pos.x=1-pos.x;
		if (pos.y<0) pos.y=1-pos.y;
		if (pos.z<0) pos.z=1-pos.z;
		return (((int)pos.x+(int)pos.y+(int)pos.z)%2==0 ? black : white)->getSurface(pos,n);
	}
};

class ScalarField {
public:
	virtual float getValue(Vector3D p)=0;
};

class Noise : public ScalarField {
public:
	virtual float getValue(Vector3D p) {
		p*=10;
		return (cos(p.x)+cos(p.y)+cos(p.z))/3;
	}
};

class Dalle : public ScalarField {
private:
	float size;
public:
	Dalle(float s) {
		size=s;
	}
	virtual float getValue(Vector3D p) {
		p/=size;
		if (p.x<0) p.x=1-p.x;
		if (p.y<0) p.y=1-p.y;
		if (p.z<0) p.z=1-p.z;
		float fx=p.x-int(p.x);
		float fy=p.y-int(p.y);
		float m=0.1;
		return min(0.05,min(fx,min(fy,min(1-fx,1-fy))));
	}
};

class HeightMapper : public Material {
private:
	Material *material;
	ScalarField *heightField;
	float scale;
public:
	HeightMapper(Material *m,ScalarField *f,float s) {
		material=m;
		heightField=f;
		scale=s;
	}
	virtual Surface getSurface(Vector3D pos,Vector3D n) {
		Surface s=material->getSurface(pos,n);
		Vector3D u=s.n^Vector3D(1,0.4,0.87); 
		if (u.length()<1e-4)
			u=s.n^Vector3D(0,-0.4,0.9); 
		u.normalize();
		Vector3D v=s.n^u; v.normalize();
		const float d=0.01;
		float dh_du=(heightField->getValue(pos+d*u)-heightField->getValue(pos-d*u))/(2*d);
		float dh_dv=(heightField->getValue(pos+d*v)-heightField->getValue(pos-d*v))/(2*d);
		s.n-=(dh_du*u+dh_dv*v)*scale;
		s.n.normalize();
		return s;
	}
};


class TransformedMaterial : public Material {
private:
	Matrix invTransform;
	Material *material;
public:
	TransformedMaterial(Matrix t,Material *m) {
		invTransform=t;
		invTransform.invert();
		material=m;
	}
	virtual Surface getSurface(Vector3D pos,Vector3D n) {
		return material->getSurface(invTransform*pos,n);
	}
};

struct Intersection {
	float dist;
	Vector3D location;
	Vector3D n;
	Material *material;
};

class Object {
public:
	bool transparent;
	Material *material;
	Intersection *pool;
	bool inside;
	virtual int intersections(const Ray &r)=0; 
	// Есть побочный эффект - устанавливает флаг inside,
	// он указывает, внутри или снаружи находится начало луча
	Object() {
		transparent=false;
	}
};

class Plane : public Object {
private:
	Vector3D n;
	float d;
public:
	Plane(Vector3D n,float d,Material *m=NULL) {
		this->n=n/n.length();
		this->d=d/n.length();
		material=m;
		pool=new Intersection[1];
	}
	virtual int intersections(const Ray &r) {
		float q=(r.origin&n)-d;
		inside=q<0;
		float k=r.dir&n;
		if (abs(k)<1e-5) return 0;
		if (k>=0 && q>=0 || k<=0 && q<=0) return 0;
		float dist=-q/k;
		pool->dist=dist;
		pool->location=r.origin+dist*r.dir;
		pool->n=n;
		pool->material=material;
		return 1;
	}
};

class Sphere : public Object {
private:
	Vector3D center;
	float r;
public:
	Sphere(Vector3D c,float r,Material *m=NULL) {
		center=c; this->r=r;
		pool=new Intersection[2];
		material=m;
	}
	virtual int intersections(const Ray &r) {
		float p=-(r.origin-center)&r.dir;
		float q=((r.origin-center)&(r.origin-center))-this->r*this->r;
		inside=q<0;
		if (p<0 && q>0)
			return 0;
		if (p*p<q+1e-6)
			return 0;
		float sqrt_d_4=sqrt(p*p-q);
		int count=0;
		if (q>0) {
			pool[count].dist=p-sqrt_d_4;
			pool[count].location=r.origin+pool[count].dist*r.dir;
			pool[count].n=(1/this->r)*(pool[count].location-center);
			pool[count].material=material;
			count++;
		}
		pool[count].dist=p+sqrt_d_4;
		pool[count].location=r.origin+pool[count].dist*r.dir;
		pool[count].n=(1/this->r)*(pool[count].location-center);
		pool[count].material=material;
		count++;
		return count;
	}
};

class LogOp : public Object {
private:
	Object *a,*b;
	bool singleMaterial;
	bool and; // для того, чтобы можно было получать как пересечение, так и объединение (по двойственности)
public:
	LogOp(Object *a,Object *b,bool and) {
		this->a=a; this->b=b;
		singleMaterial=false;
		this->and=and;
		pool=new Intersection[100];
	}
	LogOp(Object *a,Object *b,bool and,Material *m) {
		this->a=a; this->b=b;
		singleMaterial=true;
		material=m;
		this->and=and;
		pool=new Intersection[100];
	}
	virtual int intersections(const Ray &r) {
		int n1=a->intersections(r);
		bool ia=a->inside;
		if (ia!=and && n1==0) {
			inside=!and;
			return 0;
		}

		int n2=b->intersections(r);
		bool ib=b->inside;

		inside=and? ia&&ib : ia||ib;

		int count=0;
		int i1=0,i2=0;
		while (i1<n1 && i2<n2) {
			if (a->pool[i1].dist<b->pool[i2].dist) {
				ia=!ia;
				if (ib==and) {
					pool[count]=a->pool[i1];
					count++;
				}
				i1++;
			} else {
				ib=!ib;
				if (ia==and) {
					pool[count]=b->pool[i2];
					count++;
				}
				i2++;
			}
		}
		if (ib==and) 
			while (i1<n1) {
				pool[count]=a->pool[i1];
				count++;
				i1++;
			}
		if (ia==and) 
			while (i2<n2) {
				pool[count]=b->pool[i2];
				count++;
			i2++;
			}
		if (singleMaterial) 
			for (int i=0;i<count;i++)
				pool[i].material=material;
		return count;
	}
};

Object *and(Object *a,Object *b,Material *m=NULL) {
	if (m==NULL) return new LogOp(a,b,true);
	return new LogOp(a,b,true,m);
}

Object *or(Object *a,Object *b,Material *m=NULL) {
	if (m==NULL) return new LogOp(a,b,false);
	return new LogOp(a,b,false,m);
}

class TurnedInsideOut : public Object {
private:
	Object *a;
public:
	TurnedInsideOut(Object *a) {
		this->a=a;
		pool=a->pool;
	}
	virtual int intersections(const Ray &r) {
		int n=a->intersections(r);
		inside=!a->inside;
		for (int i=0;i<n;i++)
			pool[i].n=-pool[i].n;
		return n;
	}
};

Object *turnInsideOut(Object *o) {
	return new TurnedInsideOut(o);
}

class BoundingVolume : public Object {
private:
	Object *complex,*bound;
public:
	BoundingVolume(Object *c,Object *b) {
		complex=c; bound=b;
		pool=complex->pool;
	}
	virtual int intersections(const Ray &r) {
		int n=bound->intersections(r);
		if (n==0 && !bound->inside) {
			inside=false;
			return 0;
		}
		n=complex->intersections(r);
		inside=complex->inside;
		return n;
	}
};


class TransformedObject : public Object {
private:
	Matrix transform;
	Matrix invTransform;
	Matrix3D normalTransform;
	Object *object;
public:
	TransformedObject(Matrix t,Object *o) {
		transform=t;
		invTransform=t; invTransform.invert();
		for (int i=0;i<3;i++)
			for (int j=0;j<3;j++)
				normalTransform[i][j]=t.x[i][j];
		normalTransform.invert(); normalTransform.transpose();
		object=o;
		pool=object->pool;
	}
	virtual int intersections(const Ray &r) {
		Ray ray;
		ray.origin=invTransform*r.origin;
		ray.dir=invTransform*r.dir-invTransform*Vector3D(0,0,0);
		ray.dir.normalize();
		int n=object->intersections(ray);
		inside=object->inside;
		for (int i=0;i<n;i++) {
			pool[i].location=transform*pool[i].location;
			pool[i].dist=(r.origin-pool[i].location).length();
			pool[i].n=normalTransform*pool[i].n;
		}
		return n;
	}
};

class Quadric : public Object {
private:
	Matrix3D a;
	Vector3D b;
	float c;
	//pT*a*p+b*p+c=0
public:
	Quadric(Matrix3D a, Vector3D b,float c,Material *m=NULL) {
		this->a=a; 
		this->a.transpose(); 
		this->a+=a;
		this->a*=0.5;
		this->b=b; this->c=c;
		material=m;
		pool=new Intersection[2];
	}
	Vector3D getNormal(Vector3D p) {
		Vector3D n=2*(a*p)+b;
		n.normalize();
		return n;
	}
	void makeIntersection(Intersection &i,float d,const Ray &r) {
		i.dist=d;
		i.location=r.origin+d*r.dir;
		i.material=material;
		i.n=getNormal(i.location);
	}
	virtual int intersections(const Ray &r) {
		//float a=r.dir&(this->a*r.dir);
		float a=
			this->a[0][0]*r.dir.x*r.dir.x+
			this->a[1][1]*r.dir.y*r.dir.y+
			this->a[2][2]*r.dir.z*r.dir.z+
			2*(
			this->a[0][1]*r.dir.x*r.dir.y+
			this->a[1][2]*r.dir.y*r.dir.z+
			this->a[0][2]*r.dir.x*r.dir.z
			);
		float b=2*(r.origin&(this->a*r.dir))+(this->b&r.dir);
		float c=(r.origin&(this->a*r.origin))+(this->b&r.origin)+this->c;
		inside=c<0;
		if (abs(a)<1e-5) {
			if (abs(b)<1e-4) return 0;
			if (b>0 && c>=0 || b<0 && c<=0) return 0;
			makeIntersection(pool[0],-c/b,r);
			return 1;
		}
		if (a<0) {
			a=-a; b=-b; c=-c;
		}
		if (b>0 && c>0) return 0;
		float d=b*b-4*a*c;
		if (d<=0) return 0;
		d=sqrt(d);
		int count=0;
		if (c>0) {
			makeIntersection(pool[count],(-b-d)/(2*a),r);
			count++;
		}
		makeIntersection(pool[count],(-b+d)/(2*a),r);
		count++;
		return count;
	}
};

vector<Object*> objects;

bool findClosestIntersection(const Ray &r,Intersection &intr) {
	intr.dist=1e10;
	for (int i=0;i<objects.size();i++) {
		if (objects[i]->intersections(r)>0) 
			if (intr.dist>objects[i]->pool->dist)
				intr=*objects[i]->pool;
	}
	if (intr.dist!=1e10) return true;
	return false;
}

struct Lighting {
	Vector3D influence;
	Vector3D direction;
	float dist;
};

class Light {
public:
	virtual Lighting getLighting(Vector3D point)=0;
};
class PointLight : public Light {
	Vector3D location;
	Vector3D color;
public:
	PointLight(Vector3D loc,Vector3D c) {
		location=loc; color=c;
	}
	virtual Lighting getLighting(Vector3D point) {
		Lighting st;
		Vector3D l=location-point;
		st.influence=color/((l&l)*4*3.14);
		st.direction=l; 
		st.dist=l.length();
		st.direction.normalize();
		return st;
	}
};

class ParallelLight : public Light {
	Vector3D direction;
	Vector3D color;
public:
	ParallelLight(Vector3D d,Vector3D c) {
		direction=d; 
		direction.normalize();
		color=c;
	}
	virtual Lighting getLighting(Vector3D point) {
		Lighting st;
		st.influence=color;
		st.direction=-direction; 
		st.dist=1e6;
		return st;
	}
};


vector<Light*> lights;

Vector3D shadowTest(const Ray &r,float maxDist) {
//	return 1;
	Vector3D transparancy=Vector3D(1,1,1);
	for (int i=0;i<objects.size();i++) {
		int numIntr=objects[i]->intersections(r);
		if (numIntr==0) continue;
		if (!objects[i]->transparent) 
			if (objects[i]->pool[0].dist<maxDist)
				return Vector3D(0,0,0);
			else continue;
		for (int j=0;j<numIntr;j++) {
			if (objects[i]->pool[j].dist>maxDist) break;
			Surface s=objects[i]->pool[j].material->getSurface(objects[i]->pool[j].location,objects[i]->pool[j].n);
			if (s.refractive.length()<1e-3) return Vector3D(0,0,0);
			transparancy*=s.refractive;
			transparancy*=0.5+0.5*abs(s.n&r.dir);
		}
	}
	return transparancy;
}

Vector3D traceRay(const Ray &r,int deep,Vector3D weight) {
	if (deep>14) return Vector3D(0.5,0.5,0.5);
	if ((weight&colorIntensity)<0.03) return Vector3D(0.5,0.5,0.5);
	Intersection intr;
	if (findClosestIntersection(r,intr)) {
		Surface surface=intr.material->getSurface(intr.location,intr.n);
		Vector3D color=surface.ambient;
		for (int i=0;i<lights.size();i++) {
			Lighting  lighting=lights[i]->getLighting(intr.location);
			float w=abs(lighting.direction & surface.n);
			Ray ray;
			ray.origin=intr.location;
			ray.dir=lighting.direction;
			ray.origin+=ray.dir*0.001;
			Vector3D shadow=shadowTest(ray,lighting.dist);
			if ((shadow&colorIntensity)>0.03) {
				color+=w*shadow*lighting.influence*surface.diffuse;
				Vector3D middle=-r.dir+lighting.direction; middle.normalize();
				float w=abs(middle&surface.n);
				if (w>1e-3)
					color+=exp(surface.phong*log(w))*shadow*lighting.influence*surface.specular;
			}
		}
		Ray ray;
		ray.origin=intr.location;
		ray.dir=r.dir-2*(surface.n&r.dir)*surface.n;
		ray.origin+=ray.dir*0.001;
		color+=surface.reflective*traceRay(ray,deep+1,weight*surface.reflective);

		if ((surface.refractive&colorIntensity)>0.03) {
			ray.origin=intr.location;
			float nr=(surface.n&r.dir)<0 ? 1/surface.nRefr : surface.nRefr;
			ray.dir=r.dir-(1-nr)*(surface.n&r.dir)*surface.n;
			ray.dir.normalize();
			ray.origin+=ray.dir*0.001;
			color+=surface.refractive*traceRay(ray,deep+1,weight*surface.refractive);
		}

		return color;
	}
	return Vector3D(0,0,0);
	//return r.dir.z*Vector3D(1,0,0);
}

Vector3D buffer[1024][1280];

void makePicture(int frame,int w, int h) {

	for (int j=0;j<=h;j++)
		for (int i=0;i<=w;i++)
			buffer[j][i]=traceRay(camera->getRay((float)i/w,(float)j/h),0,Vector3D(1,1,1));;

	FILE *f=_popen("\"c:\\Program Files\\python24\\python.exe\" text2bmp.py","wt");
	fprintf(f,"out%.4d.bmp\n",frame);
	fprintf(f,"%d %d\n",w,h);
	for (int j=0;j<h;j++) {
		for (int i=0;i<w;i++) {
			Vector3D c=buffer[j][i]+buffer[j][i+1]+buffer[j+1][i]+buffer[j+1][i+1]; c/=4;
			if (AA) {
				float disp=0;
				for (int qx=0;qx<1;qx++)
					for (int qy=0;qy<1;qy++) 
						disp+=((c-buffer[j+qy][i+qx])*colorIntensity)&((c-buffer[j+qy][i+qx])*colorIntensity);
				if (disp>0.004) {
					//c=Vector3D(0,100,0);
					for (int k=4;k<30;k++) {
						float qx=rand()%1000*0.001;
						float qy=rand()%1000*0.001;
						c=(1/(k+1.0))*(c*k+traceRay(camera->getRay((float)(i+qx)/w,(float)(j+qy)/h),0,Vector3D(1,1,1)));
					}
				}
			} 
			if (c.x<0) c.x=0;
			if (c.y<0) c.y=0;
			if (c.z<0) c.z=0;
			if (c.x>1) c.x=1;
			if (c.y>1) c.y=1;
			if (c.z>1) c.z=1;
			c*=255;
			fprintf(f,"%d %d %d ",(int)c.x,(int)c.y,(int)c.z);
		}
		fprintf(f,"\n");
	}
	_pclose(f);
}

Object *jaggedSphere(Material *m) {
	Object *o=new Sphere(Vector3D(0,0,0),1,m);
	for (int i=0;i<8;i++) {
		Vector3D p=Vector3D(i%2,(i/2)%2,i/4)-Vector3D(0.5,0.5,0.5);
		o=and(o,turnInsideOut(new Sphere(p*0.55,0.6,m)));
	}
	return o;
}

Object *rotative(Vector3D point,Vector3D dir,float r,float angle,Material *m=NULL) {
	dir.normalize();
	Matrix3D vecStack;
	for (int i=0;i<3;i++) {
		vecStack[i][0]=dir.x;
		vecStack[i][1]=dir.y;
		vecStack[i][2]=dir.z;
	}
	Matrix3D a=Matrix3D(1)-Matrix3D::scale(dir)*vecStack;
	Matrix3D aT=a; aT.transpose();
	a=aT*a;
	Matrix3D along=(1/sqrt(3.0))*vecStack;
	Matrix3D alongT=along; alongT.transpose();
	a-=sin(angle)/cos(angle)*alongT*along;
	Vector3D b=-2*(a*point);
	float c=(point&(a*point))-r*r;
	return new Quadric(a,b,c,m);
}

Object *cylinder(Vector3D begin,Vector3D end,float r,Material *m=NULL) {
	//return rotative(begin,end-begin,r,0,m);
/*	return 
		and(
			rotative(begin,end-begin,r,0),
			new Plane(end-begin,end&(end-begin)),
		m);*/
	return 
		and(
			rotative(begin,end-begin,r,0),
			and(
				new Plane(end-begin,end&(end-begin)),
				new Plane(begin-end,begin&(begin-end))
			),
			m);
}


Object *facetedGlass(float h1,float h2,float r1,float r2,Material *glass,float fill=0,Material *drink=NULL) {
	Object *outer=new Plane(Vector3D(1,0,0),r2);
	const int n=10;
	for (int i=1;i<n;i++) 
		outer=and(new Plane(Vector3D(cos(i*2*3.14159/n),sin(i*2*3.14159/n),0),r2),outer);

	//outer=rotative(Vector3D(0,0,0),Vector3D(0,0,1),r2,0);
	
	outer=and(outer,and(new Plane(Vector3D(0,0,1),h2),new Plane(Vector3D(0,0,-1),0)));
	Object *inner=and(
		rotative(Vector3D(0,0,0),Vector3D(0,0,1),r1,0),
		new Plane(Vector3D(0,0,-1),-h1));

	//return  and(outer,turnInsideOut(inner),m);
	Object* res=and(outer,turnInsideOut(inner),glass);
	if (fill>0)
		res=or(res,cylinder(Vector3D(0,0,h1+0.001),Vector3D(0,0,h1+0.001+fill*(h2-h1-0.001)),r1-0.001,drink));
	return  new BoundingVolume(
		res,
		new Sphere(Vector3D(0,0,h2/2),sqrt(h2*h2/4+1.2*r2*r2)));
}

Object *retort(float r,float r2,float h,float wall,Material *m=NULL) {
	Object *outside=or(new Sphere(Vector3D(0,0,0),r),cylinder(Vector3D(0,0,0),Vector3D(0,0,r+h),r2));

	Object *inside=or(new Sphere(Vector3D(0,0,0),r-wall),cylinder(Vector3D(0,0,0),Vector3D(0,0,r+h+wall),r2-wall));
	return and(outside,turnInsideOut(inside),m);
	//return outside;
}

Object *testTube(float r,float h,float wall,Material *glass,float fill=0,Material *drink=NULL) {
	Object *outside=
		or(
			cylinder(Vector3D(0,0,r),Vector3D(0,0,h),r),
			new Sphere(Vector3D(0,0,r),r));
	Object *inside=
		or(
		cylinder(Vector3D(0,0,r),Vector3D(0,0,h+wall),r-wall),
		new Sphere(Vector3D(0,0,r),r-wall));
	Object *res=and(outside,turnInsideOut(inside),glass);
	if (fill>0) {
		Object *drinkInside=
			or(
			cylinder(Vector3D(0,0,r),Vector3D(0,0,h+wall),r-wall-0.001),
			new Sphere(Vector3D(0,0,r),r-wall-0.001));
		res=or(res,and(drinkInside,new Plane(Vector3D(0,0,1),h*fill),drink));
	}
	return res;
}

Object *box(Vector3D p1,Vector3D p2,Material *m=NULL) {
	return new BoundingVolume(
		and(
			new Plane(Vector3D(1,0,0),p2.x),
			and(
				new Plane(Vector3D(-1,0,0),-p1.x),
				and(
					new Plane(Vector3D(0,1,0),p2.y),
					and(
						new Plane(Vector3D(0,-1,0),-p1.y),
						and(
							new Plane(Vector3D(0,0,1),p2.z),
							new Plane(Vector3D(0,0,-1),-p1.z)
						)
					)
				)
			),
			m
		),
		new Sphere(0.5*(p1+p2),0.5*(p2-p1).length()));
}

Material *redPlastic=new UniformMaterial(Black,Red,Black,0,0.09*White);
Material *bluePlastic=new UniformMaterial(Black,Blue,Black,0,0.09*White);
Material *yellowPlastic=new UniformMaterial(Black,0.7*Yellow,Yellow,50,0.09*White);


Material *paper=new UniformMaterial(Black,White);
Material *mirror=new UniformMaterial(Black,Black,White,50,White);
Material *marble=new UniformMaterial(Black,White,White,50,0.5*White);
Material *dalle=new HeightMapper(marble,new Dalle(1),1);




Material *glass=new UniformMaterial(Black,0.01*White,White,200,0.1*White,0.99*White,1.6);
Material *plumGlass=new UniformMaterial(Black,Black,2*White,200,0.1*White,Plum,1.6);
Material *greenDrink=new UniformMaterial(Black,0.1*Green,1.5*Green,60,0.1*White,0.99*Green,1.3);
Material *water=new UniformMaterial(Black,0.1*White,1.5*White,60,0.1*White,0.6*Vector3D(0.9,0.9,1),1.35);

Material *redDrink=new UniformMaterial(Black,0.1*Red,1*Red,100,0.1*White,0.99*Red,1.3);


Material *metal=new UniformMaterial(Black,0.3*White,White,200,0.6*White,Black,0);

Object *tubeBattery() {
	Object *res=
		or(
		box(Vector3D(-0.55,-0.15,0),Vector3D(0.55,0.15,0.03)),
			or(
				box(Vector3D(-0.55,-0.15,1),Vector3D(0.55,-0.12,1.03)),
				box(Vector3D(-0.55,0.12,1),Vector3D(0.55,0.15,1.03))
				),
		redPlastic);


	//testTube(0.12,1.3,0.008,glass);

	res=or (res,
		new TransformedObject(
		Matrix::translate(Vector3D(-0.3,0,0.03)),
		testTube(0.12,1.3,0.008,glass)
		)
		);

	res=or (res,
		new TransformedObject(
		Matrix::translate(Vector3D(0,0,0.03)),
		testTube(0.12,1.3,0.008,glass,0.7,redDrink)
		)
		);

	res=or (res,
		new TransformedObject(
		Matrix::translate(Vector3D(0.3,0,0.03)),
		testTube(0.12,1.3,0.008,glass,0.5,water)
		)
		);
	res=new BoundingVolume(res,
			new Sphere(Vector3D(0,0,0.6),1.2));
	return res;

}


int main(int argc, char* argv[])
{
	const Vector3D period=Vector3D(18,-2,0);
	int totalFrames=atoi(argv[2]);
	if (totalFrames==0) totalFrames=1;
	int frame=atoi(argv[1]);

//	totalFrames=100;
	//frame=13;

	float progress=(float)frame/totalFrames;

	objects.push_back(new Plane(Vector3D(0,0,1),0,dalle));

	objects.push_back(and(
		rotative(Vector3D(-1.5,-1,0.5),Vector3D(0,0,1),0,0.4),
		new Plane(Vector3D(0,0,1),0.25),
		yellowPlastic
		));

	objects.push_back(
		facetedGlass(0.05,1.3,0.47,0.5,glass,0.5,greenDrink)
	);
	objects.back()->transparent=true;
	

	objects.push_back(
		new TransformedObject(
			Matrix::translate(Vector3D(3,-1.3,1.002)),
			retort(1,0.3,0.3,0.04,glass))
		);
	objects.back()->transparent=true;

	objects.push_back(new TransformedObject(
		Matrix::translate(Vector3D(5.6,-0.8,0))*Matrix::rotateZ(2.0)*Matrix::translate(Vector3D(0,0,0.121))*Matrix::rotateX(3.14159/2)*Matrix::translate(Vector3D(0,0,-0.12)),
		testTube(0.12,1.3,0.008,glass)));
	objects.back()->transparent=true;

	objects.push_back(new TransformedObject(
		Matrix::translate(Vector3D(8.3,-1.5,0))*Matrix::rotateZ(0.7)*Matrix::translate(Vector3D(0,0,0.121))*Matrix::rotateX(3.14159/2)*Matrix::translate(Vector3D(0,0,-0.12)),
		testTube(0.12,1.3,0.008,glass)));
	objects.back()->transparent=true;


	objects.push_back(new TransformedObject(
		Matrix::translate(Vector3D(10,-1.5,0))*Matrix::rotateZ(0.5),
		tubeBattery())
	);
	objects.back()->transparent=true;


	for (int i=0;i<2;i++) {
		objects.push_back(cylinder(Vector3D(-3,-0.2,0.1)+period*i,Vector3D(-2.2,0.4,0.1)+period*i,0.1,bluePlastic));
		objects.push_back(cylinder(Vector3D(-5,0.2,0.1)+period*i,Vector3D(-4.0,-0.2,0.1)+period*i,0.1,redPlastic));
		objects.push_back(cylinder(Vector3D(-5.1,-0.7,0.1)+period*i,Vector3D(-4.6,0.2,0.3)+period*i,0.1,metal));
	}

	lights.push_back(new PointLight(Vector3D(1,-3,8),(1-progress)*700*White));
	lights.push_back(new PointLight(Vector3D(1,-3,8)+period,progress*700*White));

//	lights.push_back(new PointLight(Vector3D(3,-3,8),700*White));
//	lights.push_back(new ParallelLight(Vector3D(-1,0.5,-1),White));
//	lights.push_back(new PointLight(Vector3D(-3,-1,2),60*White));

	camera=new Camera(
		Vector3D(-5.03,-2.52,3.401)+period*progress,
		Vector3D(-4.31,0,0.51)+period*progress,
		1,0.75,Vector3D(0.01,0,1));

	float beginTime=GetTickCount();
	makePicture(frame,640,480);
	float time=0.001*(GetTickCount()-beginTime);
	printf("%0.1f\n",time);
	//scanf("%*s");
	return 0;
}

