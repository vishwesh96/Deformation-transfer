#include "Mesh.hpp"
#include <cstdlib>
#include <vector>
#include "DGP/Plane3.hpp"
#include <iostream>
#include <fstream>

using namespace std;

void calcuateS(Mesh & source_mesh, Mesh & deformed_source_mesh, vector<Matrix3> & S);
void calcuateA_c(vector<pair<long,MeshFace *> > &M,vector<Matrix3> &S, Mesh & target_mesh);
void calculateM(string correspondence_path, Mesh & source_mesh, Mesh & target_mesh, vector<pair<long,MeshFace *> > &M);


int
usage(int argc, char * argv[])
{
  DGP_CONSOLE << "";
  DGP_CONSOLE << "Usage: " << argv[0] << " <source-in> <source-deformed> <target-in> <source_target_correspondence>";
  DGP_CONSOLE << "";

  return -1;
}

int
main(int argc, char * argv[])
{
  if (argc < 2)
    return usage(argc, argv);

  if (argc == 3){
    std::string path = argv[1];
    std::string out_path = argv[2];
 
   Mesh mesh ;
    if (!mesh.load(path)){
      return -1;
    }
    std::ofstream out(out_path.c_str(), std::ios::binary);
    if (!out)
    {
      DGP_ERROR << "Could not open '" << out_path << "' for writing";
      return -1;
    }

   for(auto it = mesh.facesBegin(); it!=mesh.facesEnd(); ++it){
      auto vit = (*it).verticesBegin(); 
      Vector3 v0 = (*vit)->getPosition();
      vit++;
      Vector3 v1 = (*vit)->getPosition();
      vit++;
      Vector3 v2 = (*vit)->getPosition(); 
      Vector3 point  = (v0+v1+v2)/3;
      Vector3 normal = (*it).getNormal();
      out << point[0] << ' ' << point[1] << ' ' << point[2] << ' ' <<normal[0] << ' ' << normal[1] << ' ' << normal[2] << '\n';
    }
  return 0;
  }

  std::string source_path = argv[1];
  std::string source_deformed_path = argv[2];
  std::string target_path = argv[3];
  std::string correspondence_path = argv[4];

  Mesh source_mesh, deformed_source_mesh, target_mesh;
  
  if (!source_mesh.load(source_path) || !deformed_source_mesh.load(source_deformed_path) || !target_mesh.load(target_path))
    return -1;

  DGP_CONSOLE << "Read mesh '" << source_mesh.getName() << "' with " << source_mesh.numVertices() << " vertices, " << source_mesh.numEdges()
              << " edges and " << source_mesh.numFaces() << " faces from " << source_path;
  DGP_CONSOLE << "Read mesh '" << deformed_source_mesh.getName() << "' with " << deformed_source_mesh.numVertices() << " vertices, " << deformed_source_mesh.numEdges()
              << " edges and " << deformed_source_mesh.numFaces() << " faces from " << source_deformed_path;
  DGP_CONSOLE << "Read mesh '" << target_mesh.getName() << "' with " << target_mesh.numVertices() << " vertices, " << target_mesh.numEdges()
              << " edges and " << target_mesh.numFaces() << " faces from " << target_path;


  //setting id of meshvertex
  long li = 0;
  for(auto it = source_mesh.verticesBegin(); it!=source_mesh.verticesEnd(); ++it, li++){it->id = li;
    if((it)->facesBegin()==(it)->facesEnd())
      DGP_CONSOLE << "\nSOURCE NO FACE" <<" "<<li;
  }

  li = 0;
  for(auto it = deformed_source_mesh.verticesBegin(); it!=deformed_source_mesh.verticesEnd(); ++it, li++)
  { it->id = li;
    if((it)->facesBegin()==(it)->facesEnd())
      DGP_CONSOLE << "\nDEFORMED SOURCE NO FACE" <<" "<<li;
  }

  li = 0;
  for(auto it = target_mesh.verticesBegin(); it!=target_mesh.verticesEnd(); ++it, li++)
  {
    it->id = li;
    if((it)->facesBegin()==(it)->facesEnd())
      DGP_CONSOLE << "\nTARGET NO FACE" <<" "<<li;
  }


  vector<Matrix3> S;								//source deformations
  vector<pair<long,MeshFace *> > M;
  calculateM(correspondence_path,source_mesh,target_mesh,M);
  // auto it = target_mesh.facesBegin();
  
  // for(int i=0;it!=target_mesh.facesEnd();i++,it++)
  // {
  //   M.push_back(pair<long,MeshFace*>(i,&(*it)));
  // }

  calcuateS(source_mesh,deformed_source_mesh,S);	//calcuate S
  calcuateA_c(M,S,target_mesh);
  return 0;
}


void calcuateS(Mesh & source_mesh, Mesh & deformed_source_mesh, vector<Matrix3> & S)
{
	for(Mesh::FaceIterator fit = source_mesh.facesBegin(), dfit = deformed_source_mesh.facesBegin(); fit!= source_mesh.facesEnd() && dfit!=deformed_source_mesh.facesEnd(); fit++, dfit++){
		list<Mesh::Vertex *>::iterator  vit = fit->verticesBegin();
		Vector3 v1 = (*vit)->getPosition();
		vit++;
		Vector3 v2 = (*vit)->getPosition();
		vit++;
		Vector3 v3 = (*vit)->getPosition();
		Vector3 v4 = v1 + ((v2 - v1).cross((v3 - v1)))/sqrt(((v2 - v1).cross(v3 - v1)).length());

		list<Mesh::Vertex *>::iterator dvit = dfit->verticesBegin();
		Vector3 dv1 = (*dvit)->getPosition();
		dvit++;
		Vector3 dv2 = (*dvit)->getPosition();
		dvit++;
		Vector3 dv3 = (*dvit)->getPosition();
		Vector3 dv4 = dv1 + ((dv2 - dv1).cross((dv3 - dv1)))/sqrt(((dv2 - dv1).cross(dv3 - dv1)).length());

		Matrix3  V = Matrix3::fromColumns(v2-v1,v3-v1,v4-v1);  
		Matrix3  dV = Matrix3::fromColumns(dv2-dv1,dv3-dv1,dv4-dv1);  

		if(V.determinant()<1e-10){
			DGP_CONSOLE<<"determinant error"<<std::endl;
			exit(0);
		}
		S.push_back(dV*V.inverse());
	}
}

void calcuateA_c(vector<pair<long,MeshFace*> > &M,vector<Matrix3> &S, Mesh & target_mesh)
{

  ofstream a_file;
  ofstream c_file;
  a_file.open("A.txt");
  c_file.open("c.txt");
  a_file << 9*M.size() << " " << 3*(target_mesh.numVertices()+target_mesh.numFaces()) << "\n";
  c_file << 9*M.size() << " " << 1 << "\n"; 
	for(unsigned long i=0; i < M.size(); i++)
	{
		list<Mesh::Vertex *>::iterator  vit = M[i].second->verticesBegin();
		Vector3 v1 = (*vit)->getPosition();
		MeshVertex *vert1 = (*vit);
		vit++;
		Vector3 v2 = (*vit)->getPosition();
		MeshVertex *vert2 = (*vit);
		vit++;
		Vector3 v3 = (*vit)->getPosition();
		MeshVertex *vert3 = (*vit);
		Vector3 v4 = v1 + ((v2 - v1).cross((v3 - v1)))/sqrt(((v2 - v1).cross(v3 - v1)).length());
		Matrix3 V = Matrix3::fromColumns(v2-v1,v3-v1,v4-v1);  

		if(V.determinant()<1e-10)
		{
			DGP_CONSOLE<<"determinant error"<<std::endl;
			exit(0);
		}
		V.invert();

		for(int j=0;j<3;j++)
		{
			for(int k=0;k<3;k++)
			{
        long index = i*9+j*3+k;
        c_file << index << " " << S[M[i].first](j,k) << "\n";
        a_file << index << " " << vert1->id*3+j << " " << -(V(0,k)+V(1,k)+V(2,k)) << "\n";
        a_file << index << " " << vert2->id*3+j << " " << V(0,k) << "\n";
        a_file << index << " " << vert3->id*3+j << " " << V(1,k) << "\n";
        a_file << index << " " << target_mesh.numVertices()*3 + M[i].second->id*3 + j << " " << V(2,k) << "\n";
			}
		}
	}
  a_file.close();
  c_file.close();
}


void calculateM(string correspondence_path, Mesh & source_mesh, Mesh & target_mesh, vector<pair<long,MeshFace *> > &M)
{
    std::ifstream in(correspondence_path.c_str());
    if (!in)
    {
      DGP_ERROR << "Could not open correspondence map" << correspondence_path << "' for reading";
      return;
    } 
  long np;
    if (!(in >>np))
    {
      DGP_ERROR << "Could not read number of correspondence points from file '" << correspondence_path << '\'';
      return;
    }
    vector<Vector3> source_vertices;
    Vector3 s;
    Vector3 t;
    for (long i = 0; i < np; ++i)
    {
      if (!(in >> t[0] >> t[1] >> t[2] >> s[0] >> s[1] >> s[2]))
      {
        DGP_ERROR << "Could not read correspondence points ";
      }
      source_vertices.push_back(s);
    }
    // if(target_mesh.numFaces() != np){
    //     DGP_ERROR<<"Number of correspondence points is not equal to number of faces of target mesh";
    //     return ;
    // }
    auto tit = target_mesh.facesBegin();
    for(long i=0;i<np;i++,tit++){
      long source_index = 0;
      MeshFace * target_face = &(*tit) ;
      double min_dist = numeric_limits<double>::max();
      for(auto it = source_mesh.facesBegin(); it!= source_mesh.facesEnd(); it++){
          auto vit = (*it).verticesBegin(); 
          Vector3 v0 = (*vit)->getPosition();
          vit++;
          Vector3 v1 = (*vit)->getPosition();
          vit++;
          Vector3 v2 = (*vit)->getPosition();
          // double temp_dist = (v0 - source_vertices[i]).length() + (v1 - source_vertices[i]).length() + (v2 - source_vertices[i]).length();
          double temp_dist = ((v0+v1+v2)/3.0 - source_vertices[i]).length();
          if(temp_dist < min_dist){
            min_dist = temp_dist;
            source_index = (*it).id;
          }
        }     
      M.push_back(pair<long,MeshFace*>(source_index,target_face));
    }   
  // auto tit = target_mesh.facesBegin();
  //   for(long i=0;i<np;i++,tit++){
  //     long source_index = 0;
  //     MeshFace * target_face = &(*tit) ;
  //     double min_dist = numeric_limits<double>::max();
  //     for(auto it = source_mesh.facesBegin(); it!= source_mesh.facesEnd(); it++){
  //         auto vit = (*it).verticesBegin(); 
  //         Vector3 v1 = (*vit)->getPosition();
  //         vit++;
  //         Vector3 v2 = (*vit)->getPosition();
  //         vit++;
  //         Vector3 v3 = (*vit)->getPosition();
  //         Vector3 proj = source_vertices[i] + ((v1-source_vertices[i]).dot((*it).getNormal()))*(*it).getNormal();
          
  //         double u = ((proj[0]*v2[1])-(proj[0]*v3[1]) - (v2[0]*proj[1]) + (v2[0]*v3[1])  + (v3[0]*proj[1]) - (v3[0]*v2[1]))/ ((v1[0] * v2[1])  - (v1[0] * v3[1])  - (v2[0] * v1[1]) + (v2[0] * v3[1]) + (v3[0] * v1[1])  - (v3[0] * v2[1]));
  //         double v = ((v1[0] * proj[1]) - (v1[0] * v3[1]) - (proj[0] * v1[1]) + (proj[0] * v3[1]) + (v3[0] * v1[1]) - (v3[0] * proj[1]))/((v1[0] * v2[1])  - (v1[0] * v3[1])  - (v2[0] * v1[1]) + (v2[0] * v3[1]) + (v3[0] * v1[1])  - (v3[0] * v2[1]));
  //         double w = ((v1[0] * v2[1]) - (v1[0] * proj[1]) - (v2[0] * v1[1]) + (v2[0] * proj[1]) + (proj[0] * v1[1]) - (proj[0] * v2[1]))/((v1[0] * v2[1])  - (v1[0] * v3[1])  - (v2[0] * v1[1]) + (v2[0] * v3[1]) + (v3[0] * v1[1])  - (v3[0] * v2[1]));

  //         Vector3 nearest =  Vector3(u,v,w);
  //         nearest.unitize();

  //         nearest = v1*nearest[0] + v2*nearest[1] + v3*nearest[2];

  //         double temp_dist = (source_vertices[i] - nearest).length();
  //         if(temp_dist < min_dist){
  //           min_dist = temp_dist;
  //           source_index = (*it).id;
  //         }
  //       }     
  //     M.push_back(pair<long,MeshFace*>(source_index,target_face));
  //   }   
}