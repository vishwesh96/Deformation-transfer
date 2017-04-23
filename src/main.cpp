#include "Mesh.hpp"
#include "Viewer.hpp"
#undef Success  
#include "Eigen/Sparse"
#include <cstdlib>
#include <vector>


using namespace std;

void calcuateS(Mesh & source_mesh, Mesh & deformed_source_mesh, vector<Matrix3> & S);
void calcuateA_c(vector<pair<long,MeshFace *> > &M,vector<Matrix3> &S, Mesh & target_mesh, Eigen::SparseMatrix<double,Eigen::RowMajor> & A, Eigen::SparseVector<double> & c);

int
usage(int argc, char * argv[])
{
  DGP_CONSOLE << "";
  DGP_CONSOLE << "Usage: " << argv[0] << " <source-in> <source-deformed> <target-in>";
  DGP_CONSOLE << "";

  return -1;
}

int
main(int argc, char * argv[])
{
  if (argc < 2)
    return usage(argc, argv);

  std::string source_path = argv[1];
  std::string source_deformed_path = argv[2];
  std::string target_path = argv[3];

  Mesh source_mesh, deformed_source_mesh, target_mesh, deformed_target_mesh;
  
  if (!source_mesh.load(source_path) || !deformed_source_mesh.load(source_deformed_path) || !target_mesh.load(target_path) || !deformed_target_mesh.load(target_path))
    return -1;

  DGP_CONSOLE << "Read mesh '" << source_mesh.getName() << "' with " << source_mesh.numVertices() << " vertices, " << source_mesh.numEdges()
              << " edges and " << source_mesh.numFaces() << " faces from " << source_path;
  DGP_CONSOLE << "Read mesh '" << deformed_source_mesh.getName() << "' with " << deformed_source_mesh.numVertices() << " vertices, " << deformed_source_mesh.numEdges()
              << " edges and " << deformed_source_mesh.numFaces() << " faces from " << source_deformed_path;
  DGP_CONSOLE << "Read mesh '" << target_mesh.getName() << "' with " << target_mesh.numVertices() << " vertices, " << target_mesh.numEdges()
              << " edges and " << target_mesh.numFaces() << " faces from " << target_path;

  // if (target_num_faces >= 0 && mesh.numFaces() > target_num_faces)
  // {
  //   mesh.updateBounds();
  // }

  // if (!out_path.empty())
  // {
  //   if (!mesh.save(out_path))
  //     return -1;

  //   DGP_CONSOLE << "Saved mesh to " << out_path;
  // }

  // source_mesh.check(deformed_source_mesh);			//check source and deformed dimensions and also check for triangular mesh

  vector<Matrix3> S;								//source deformations
  vector<pair<long,MeshFace *> > M;
  // M.resize(10);
  auto it = target_mesh.facesBegin();
  
  DGP_CONSOLE<<"hi\n";
  for(int i=0;i<10;i++,it++)
  {
  	DGP_CONSOLE<<i<<endl;
  	M.push_back(pair<long,MeshFace*>(i,&(*it)));
  }
  DGP_CONSOLE<<"bye\n";

  Eigen::SparseMatrix<double,Eigen::RowMajor> A(9*M.size(),3*(target_mesh.numVertices()+target_mesh.numFaces())); 
  Eigen::SparseVector<double> c(9*M.size());
  calcuateS(source_mesh,deformed_source_mesh,S);	//calcuate S
  calcuateA_c(M,S,target_mesh,A,c);

  //DGP_CONSOLE<<A;

  Viewer viewer1, viewer2, viewer3;
  viewer1.setObject(&source_mesh);
  viewer1.launch(argc, argv);

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
		Vector3 v4 = v1 + ((v2 - v1).cross((v3 - v1)))/((v2 - v1).cross(v3 - v1)).length();

		list<Mesh::Vertex *>::iterator dvit = dfit->verticesBegin();
		Vector3 dv1 = (*dvit)->getPosition();
		dvit++;
		Vector3 dv2 = (*dvit)->getPosition();
		dvit++;
		Vector3 dv3 = (*dvit)->getPosition();
		Vector3 dv4 = dv1 + ((dv2 - dv1).cross((dv3 - dv1)))/((dv2 - dv1).cross(dv3 - dv1)).length();

		Matrix3  V = Matrix3::fromColumns(v2-v1,v3-v1,v4-v1);  
		Matrix3  dV = Matrix3::fromColumns(dv2-dv1,dv3-dv1,dv4-dv1);  

		if(V.determinant()<1e-10){
			DGP_CONSOLE<<"determinant error"<<std::endl;
			exit(0);
		}
		S.push_back(dV*V.inverse());
	}
}

void calcuateA_c(vector<pair<long,MeshFace*> > &M,vector<Matrix3> &S, Mesh & target_mesh, Eigen::SparseMatrix<double,Eigen::RowMajor> & A, Eigen::SparseVector<double> & c)
{
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
		Vector3 v4 = v1 + ((v2 - v1).cross((v3 - v1)))/((v2 - v1).cross(v3 - v1)).length();
		Matrix3 V = Matrix3::fromColumns(v2-v1,v3-v1,v4-v1);  

		if(V.determinant()<1e-10)
		{
			DGP_CONSOLE<<"determinant error"<<std::endl;
			exit(0);
		}

		V.invert();

		A.reserve(Eigen::VectorXi::Constant(A.rows(),12));
		for(int j=0;j<3;j++)
		{
			for(int k=0;k<3;k++)
			{
				// cout<<i<<" "<<j<<" "<<k<<endl;
				// [v1 v2 v3] * for all vertices then v4 corresponding to all triangles
				c.coeffRef(i*9+j*3+k) = S[M[i].first](j,k);
				A.insert(i*9+j*3+k,vert1->id*3+j) = -(V(0,k)+V(1,k)+V(2,k));
				A.insert(i*9+j*3+k,vert2->id*3+j) = V(0,k);
				A.insert(i*9+j*3+k,vert3->id*3+j) = V(1,k);
				A.insert(i*9+j*3+k,target_mesh.numVertices()*3 + M[i].second->id*3 + j) = V(2,k);
			}
		}


	}
}