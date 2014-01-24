/////////////////////////////////////////////////////////
// Sample program to book                              //
//  Computer Graphics : Dynamics & Realistic Imaging.  //
//      by A.V. Boreskoff, E.V. Shikin                 //
//                                                     //
// Author:                                             //
//    Alex V. Boreskoff                                //
//                                                     //
// E-mail:                                             //
//    alex@garser.msk.su                               //
/////////////////////////////////////////////////////////

#ifndef	__COLORS__
#define	__COLORS__

#ifndef	__VECTOR__
#include	"vector3d.h"
#endif

#define	Aquamarine		Vector3D ( 0.439216, 0.858824, 0.576471 )
#define	Black           	Vector3D ( 0,        0,        0        )
#define	Blue            	Vector3D ( 0,        0,        1        )
#define	BlueViolet      	Vector3D ( 0.623529, 0.372549, 0.623529 )
#define	Brown           	Vector3D ( 0.647059, 0.164706, 0.164706 )
#define	CadetBlue       	Vector3D ( 0.372549, 0.623529, 0.623529 )
#define	Coral           	Vector3D ( 1,        0.498039, 0        )
#define	CornflowerBlue  	Vector3D ( 0.258824, 0.258824, 0.435294 )
#define	Cyan            	Vector3D ( 0,        1,        1        )
#define	DarkGreen       	Vector3D ( 0.184314, 0.309804, 0.184314 )
#define	DarkOliveGreen  	Vector3D ( 0.309804, 0.309804, 0.184314 )
#define	DarkOrchid      	Vector3D ( 0.6,      0.196078, 0.8      )
#define	DarkSlateBlue   	Vector3D ( 0.419608, 0.137255, 0.556863 )
#define	DarkSlateGray   	Vector3D ( 0.184314, 0.309804, 0.309804 )
#define	DarkSlateGrey   	Vector3D ( 0.184314, 0.309804, 0.309804 )
#define	DarkTurquoise   	Vector3D ( 0.439216, 0.576471, 0.858824 )
#define	DimGray         	Vector3D ( 0.329412, 0.329412, 0.329412 )
#define	DimGrey         	Vector3D ( 0.329412, 0.329412, 0.329412 )
#define	Firebrick       	Vector3D ( 0.9, 0.4, 0.3 )	//Vector3D ( 0.556863, 0.137255, 0.137255 )
#define	ForestGreen     	Vector3D ( 0.137255, 0.556863, 0.137255 )
#define	Gold            	Vector3D ( 0.8,      0.498039, 0.196078 )
#define	Goldenrod       	Vector3D ( 0.858824, 0.858824, 0.439216 )
#define	Gray            	Vector3D ( 0.752941, 0.752941, 0.752941 )
#define	Green           	Vector3D ( 0,        1,        0        )
#define	GreenYellow     	Vector3D ( 0.576471, 0.858824, 0.439216 )
#define	Grey            	Vector3D ( 0.752941, 0.752941, 0.752941 )
#define	IndianRed       	Vector3D ( 0.309804, 0.184314, 0.184314 )
#define	Khaki           	Vector3D ( 0.623529, 0.623529, 0.372549 )
#define	LightBlue       	Vector3D ( 0.74902,  0.847059, 0.847059 )
#define	LightGray       	Vector3D ( 0.658824, 0.658824, 0.658824 )
#define	LightGrey       	Vector3D ( 0.658824, 0.658824, 0.658824 )
#define	LightSteelBlue  	Vector3D ( 0.560784, 0.560784, 0.737255 )
#define	LimeGreen       	Vector3D ( 0.196078, 0.8,      0.196078 )
#define	Magenta         	Vector3D ( 1,        0,        1        )
#define	Maroon          	Vector3D ( 0.556863, 0.137255, 0.419608 )
#define	MediumAquamarine 	Vector3D ( 0.196078, 0.8,      0.6      )
#define	MediumBlue      	Vector3D ( 0.196078, 0.196078, 0.8      )
#define	MediumForestGreen 	Vector3D ( 0.419608, 0.556863, 0.137255 )
#define	MediumGoldenrod   	Vector3D ( 0.917647, 0.917647, 0.678431 )
#define	MediumOrchid      	Vector3D ( 0.576471, 0.439216, 0.858824 )
#define	MediumSeaGreen    	Vector3D ( 0.258824, 0.435294, 0.258824 )
#define	MediumSlateBlue   	Vector3D ( 0.498039, 0,        1        )
#define	MediumSpringGreen 	Vector3D ( 0.498039, 1,        0        )
#define	MediumTurquoise   	Vector3D ( 0.439216, 0.858824, 0.858824 )
#define	MediumVioletRed   	Vector3D ( 0.858824, 0.439216, 0.576471 )
#define	MidnightBlue      	Vector3D ( 0.184314, 0.184314, 0.309804 )
#define	Navy              	Vector3D ( 0.137255, 0.137255, 0.556863 )
#define	NavyBlue          	Vector3D ( 0.137255, 0.137255, 0.556863 )
#define	Orange            	Vector3D ( 0.8,      0.196078, 0.196078 )
#define	OrangeRed         	Vector3D ( 0,	   0,        0.498039 )
#define	Orchid            	Vector3D ( 0.858824, 0.439216, 0.858824 )
#define	PaleGreen         	Vector3D ( 0.560784, 0.737255, 0.560784 )
#define	Pink              	Vector3D ( 0.737255, 0.560784, 0.560784 )
#define	Plum              	Vector3D ( 0.917647, 0.678431, 0.917647 )
#define	Red               	Vector3D ( 1,        0,        0        )
#define	Salmon                  Vector3D ( 0.435294, 0.258824, 0.258824 )
#define	SeaGreen                Vector3D ( 0.137255, 0.556863, 0.419608 )
#define	Sienna                  Vector3D ( 0.556863, 0.419608, 0.137255 )
#define	SkyBlue                 Vector3D ( 0.196078, 0.6,      0.8      )
#define	SlateBlue               Vector3D ( 0,        0.498039, 1        )
#define	SpringGreen             Vector3D ( 0,        1,        0.498039 )
#define	SteelBlue               Vector3D ( 0.137255, 0.419608, 0.556863 )
#define	Tan                     Vector3D ( 0.858824, 0.576471, 0.439216 )
#define	Thistle                 Vector3D ( 0.847059, 0.74902,  0.847059 )
#define	Turquoise               Vector3D ( 0.678431, 0.917647, 0.917647 )
#define	Violet                  Vector3D ( 0.309804, 0.184314, 0.309804 )
#define	VioletRed               Vector3D ( 0.8,      0.196078, 0.6      )
#define	Wheat                   Vector3D ( 0.847059, 0.847059, 0.74902  )
#define	White                   Vector3D ( 0.988235, 0.988235, 0.988235 )
#define	Yellow                  Vector3D ( 1,        1,        0        )
#define	YellowGreen             Vector3D ( 0.6,      0.8,      0.196078 )
#define	LightWood		Vector3D ( 0.6,      0.24,     0.1      )
#define	MedianWood		Vector3D ( 0.3,      0.12,     0.03     )
#define	DarkWood		Vector3D ( 0.05,     0.01,     0.005    )
#endif
