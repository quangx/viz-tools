#include <deal.II/grid/tria.h>
// Here are some functions to generate standard grids:
#include <deal.II/grid/grid_generator.h>
// Output of grids in various graphics formats:
#include <deal.II/grid/grid_out.h>

// This is needed for C++ output:
#include <iostream>
#include <fstream>
#include <cmath>
#include <random>

using namespace dealii;
enum class DataInterpretation{
  component_is_scalar,component_is_vector
};
struct StructuredData{
  private: 
     Table<3,double> priorities; 
  public:
    Table<4,double> data;
    Point<3,double> min, max;
    std::array<double,3> spacing;  
    std::array<unsigned int,3> num_values;
    bool spherical;


    StructuredData(const Point<3,double> &min,const Point<3,double> &max,const std::array<unsigned int,3> &num_values,const int &num_components,bool spherical){

        this->min=min;
        this->max=max;
        this->num_values=num_values;
        this->spherical=spherical;
        for(int i=0;i<3;++i){
          spacing[i]=(1.0*max[i]-1.0*min[i])/(1.0*num_values[i]-1);
        }
       
        TableIndices<4> t_ind(num_values[0],num_values[1],num_values[2],num_components);
        data.reinit(t_ind);
        priorities.reinit(num_values[0],num_values[1],num_values[2]);

        
    }
    StructuredData(){

    }
    Point<3,double> spherical_to_cartesian_coordinates(const std::array<double,3> &spherical_coord){
      Point<3,double> cartesian_coord;
      cartesian_coord(0)=spherical_coord[0]*std::sin(spherical_coord[2])*std::cos(spherical_coord[1]);
      cartesian_coord(1)=spherical_coord[0]*std::sin(spherical_coord[2])*std::sin(spherical_coord[1]);
      cartesian_coord(2)=spherical_coord[0]*std::cos(spherical_coord[2]);
      return cartesian_coord;

    }
    Point<3,double> spherical_to_cartesian_coordinates(Point<3,double> &spherical_coord){
      Point<3,double> cartesian_coord;
      cartesian_coord(0)=spherical_coord(0)*std::sin(spherical_coord(2))*std::cos(spherical_coord(1));
      cartesian_coord(1)=spherical_coord(0)*std::sin(spherical_coord(2))*std::sin(spherical_coord(1));
      cartesian_coord(2)=spherical_coord(0)*std::cos(spherical_coord(2));
      return cartesian_coord;

    }
    std::array<unsigned int,3> location_to_index(const Point<3,double> &p){
      std::array<unsigned int,3> index;
      for(int i=0;i<3;++i){
        unsigned int temp_index=std::floor((p(i)-min[i])/spacing[i]);
        if(((1.0*temp_index*spacing[i]+min[i])+(1.0*(temp_index+1)*spacing[i]+min[i]))/2 <= p(i)){
          ++temp_index;
        }
        index[i]=temp_index;
      }
      
      return index;
    }
    std::array<unsigned int,3> location_to_index_spherical(const Point<3,double> &p){
      std::array<unsigned int,3> index;
      for(int i=0;i<3;++i){
        unsigned int temp_index=std::floor((p(i)-min[i])/spacing[i]);
        if(((1.0*temp_index*spacing[i]+min[i])+(1.0*(temp_index+1)*spacing[i]+min[i]))/2 <= p(i)){
          ++temp_index;
        }
        index[i]=temp_index;
      }
      
      return index;
      
    }
    std::array<double,3>
    cartesian_to_spherical_coordinates(const Point<3> &position)
    {
      std::array<double,3> spherical_coord;

      spherical_coord[0] = position.norm(); // R
      spherical_coord[1] = std::atan2(position(1),position(0)); // Phi
      if (spherical_coord[1] < 0.0)
        spherical_coord[1] += 2.0*numbers::PI; // correct phi to [0,2*pi]
      
      if (spherical_coord[0] > std::numeric_limits<double>::min())
        spherical_coord[2] = std::acos(position(2)/spherical_coord[0]);
      else
        spherical_coord[2] = 0.0;
        
      return spherical_coord;
    }

    Point<3,double> index_to_location(const
    std::array<unsigned int,3> &idx){
      Point<3,double> location;
      for(int i=0;i<3;++i){
        location(i)=min[i]+idx[i]*spacing[i];
      }
      return location;
    }
    Point<3,double> index_to_location_spherical(const         //returns cartesian coordinates given an index
    std::array<unsigned int,3> &idx){
      std::array<double,3> location;
      for(int i=0;i<3;++i){
        location[i]=min[i]+idx[i]*spacing[i];
      }
      
      return spherical_to_cartesian_coordinates(location);
    }
    
    void set_values(const std::array<unsigned int,3> &idx,             
    double distance, std::vector<double> &values){
      if(1./(1e-20+distance)> priorities[idx[0]][idx[1]][idx[2]]){
        for(unsigned int i=0;i<data.size()[3];++i){  //DEBUG SHOULD BE LESS THAN VALUES.SIZE()
          
          data[idx[0]][idx[1]][idx[2]][i]=values[i];
        }
        priorities[idx[0]][idx[1]][idx[2]]=1./(1e-20+distance);
      }
    }
    std::array<unsigned int,3> approximate_extent(const Point<3,double> &position,
    double radius){
      //implement later
      std::array<unsigned int,3> result;
      for(unsigned int d=0;d<3;++d){
        result[d]=5;
      }
      return result;
    }
    
    //splat is called in write_to_vertex. in the case that spherical flag is true, p has already been converted to spherical coordinates by write_to_vertex_spherical
    void splat(const Point<3,double> &p,std::vector<double> &values,const double radius){   
      const std::array<unsigned int,3> idx=location_to_index(p); 
      std::array<unsigned int,3> extent=approximate_extent(p,radius);
      for(unsigned int iz=0;iz<2*extent[2];++iz){
        if(idx[2]+iz<extent[2]){
          continue;
        }
        for(unsigned int iy=0;iy<2*extent[1];++iy){
          if(idx[1]+iy<extent[1]){
            continue;
          }
          for(unsigned int ix=0;ix<2*extent[0];++ix){
            if(idx[0]+ix<extent[0]){
              continue;
            }
            std::array<unsigned int,3> current_index;
            current_index[0] = std::max(0l,static_cast<long>(idx[0] + ix) - static_cast<long>(extent[0]));
            current_index[1] = std::max(0l,static_cast<long>(idx[1] + iy) - static_cast<long>(extent[1]));
            current_index[2] = std::max(0l,static_cast<long>(idx[2] + iz) - static_cast<long>(extent[2]));
            if(current_index[0]>=num_values[0]||
            current_index[1]>=num_values[1] ||
            current_index[2]>=num_values[2] 
            ){
              continue;
            }
            else{
              if(!spherical){
                const double distance=index_to_location(current_index).distance(p);
                set_values(current_index,distance,values);
              }
              else{
                std::array<double,3> p_to_array={p[0],p[1],p[2]};
                
                const double distance=index_to_location_spherical(current_index).distance(spherical_to_cartesian_coordinates(p_to_array));
                set_values(current_index,distance,values);
              }
            }
          
        }
      }

    }
  }
    void to_vtk(Table<4,double> &t,const Point<3,double> min,
    const Point<3,double> max,const std::string filename,
    std::vector<DataInterpretation> component_type,const std::vector<std::string> names){
        int n1;
        int n2;
        int n3;
        std::vector<double> spacing;
        TableIndices<4> table_dim=t.size();
        
        n1=table_dim[0];
        n2=table_dim[1];
        n3=table_dim[2];
        
        for(int i=0;i<3;++i){
            spacing.push_back((double)(max[i]-min[i])/(double)(table_dim[i]-1));
        }
        
        
        
        std::ofstream file;
        file.open(filename);
        file.precision(8);
        file<<std::scientific;

        file<<"<VTKFile type=\"ImageData\" "
        "byte_order=\"LittleEndian\" header_type=\"UInt64\"> \n \t"
        "<ImageData WholeExtent=\""+std::to_string(0)+" "+std::to_string(n1-1)+""
        " "+std::to_string(0)+" "+std::to_string(n2-1)+" "+std::to_string(0)+" "+std::to_string(n3-1)+""
        "\" Origin=\"0 0 0\" Spacing=\""+std::to_string(spacing[0])+" "+std::to_string(spacing[1])+" "+std::to_string(spacing[2])+"\" Direction=\"1 0 0 0 1 0 0 0 1\">"
        "\n\t <Piece Extent=\""+std::to_string(0)+" "+std::to_string(n1-1)+" "
        ""+std::to_string(0)+" "+std::to_string(n2-1)+" "+std::to_string(0)+" "+std::to_string(n3-1)+"\">"
        "\n\t<PointData Scalars=\"T\">";
        for(unsigned int i=0;i<component_type.size();++i){
            file<<"\n\t <DataArray type=\"Float32\" Name=\""+names[i]+"\" ";
            if(component_type[i]==DataInterpretation::component_is_vector){
            file<<" NumberOfComponents=\"3\" ";

            }
            file<<"format=\"ascii\">\n";
            if(component_type[i]==DataInterpretation::component_is_vector){
            for(int z=0;z<n3;++z){
                file<<"\t";
                for(int y=0;y<n2;++y){
                for(int x=0;x<n1;++x){
                    for(int j=0;j<3;++j){
                    file<<t[x][y][z][i+j]<<" ";
                    }
                    file<<"\n\t";
                }
                file<<"\n\t";
                }
                file<<"\n";
            }
            i=i+2;
            }
            else{
            for(int z=0;z<n3;++z){
                file<<"\t";
                for(int y=0;y<n2;++y){
                for(int x=0;x<n1;++x){
                    file<<t[x][y][z][i]<<" ";
                }
                file<<"\n\t";
                }
                file<<"\n";
            }
            }
        
        
        
            file<<"</DataArray> \n ";
        }

        file<<"</PointData> \n "
        "<CellData> </CellData> \n"
        "</Piece> \n </ImageData> \n </VTKFile>";


    }
    


};

  
