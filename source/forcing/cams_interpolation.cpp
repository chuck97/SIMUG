#include "cams_interpolation.hpp"

using namespace INMOST;

namespace SIMUG
{
void CamsVectorInterpolation(SIMUG::IceMesh* mesh,
                             const NcFileInfo& file_info,
                             int netcdf_index,
                             const std::string& nc_variable_u,
                             const std::string& nc_variable_v,
                             double max_abs_value,
                             double invalid_value_fill,
                             double no_extrapolation_fill,
                             INMOST::Tag vector_tag,
                             INMOST::Tag geo_coords_tag)
{
    // get grid type
    mesh::gridType grid_type = mesh->GetMeshInfo().grid_type;

    // initialize tag with no extrapolation value
    if (grid_type == mesh::gridType::Agrid)
    {
        for(auto nodeit = mesh->GetMesh()->BeginNode(); nodeit != mesh->GetMesh()->EndNode(); ++nodeit)
        {
            nodeit->RealArray(vector_tag)[0] = no_extrapolation_fill;
            nodeit->RealArray(vector_tag)[1] = no_extrapolation_fill;
        } 
    }
    else if (grid_type == mesh::gridType::Cgrid)
    {
        for(auto edgeit = mesh->GetMesh()->BeginFace(); edgeit != mesh->GetMesh()->EndFace(); ++edgeit)
        {
            edgeit->RealArray(vector_tag)[0] = no_extrapolation_fill;
            edgeit->RealArray(vector_tag)[1] = no_extrapolation_fill;
        } 
    }
    else
    {
        SIMUG_ERR("Currently A and C grid is available for file interpolation");
    }
    BARRIER

    //Find extremal coords
    double min_geo_x;
    double min_geo_y;
    double max_geo_x;
    double max_geo_y;

    double tmp_min_x;
    double tmp_max_x;
    double tmp_min_y;
    double tmp_max_y;

    if (grid_type == mesh::gridType::Agrid)
    {
        tmp_min_x = mesh->GetMesh()->BeginNode()->RealArray(geo_coords_tag)[0]*180.0/M_PI;
        tmp_max_x = mesh->GetMesh()->BeginNode()->RealArray(geo_coords_tag)[0]*180.0/M_PI;
        tmp_min_y = mesh->GetMesh()->BeginNode()->RealArray(geo_coords_tag)[1]*180.0/M_PI;
        tmp_max_y = mesh->GetMesh()->BeginNode()->RealArray(geo_coords_tag)[1]*180.0/M_PI;

        for(auto nodeit = mesh->GetMesh()->BeginNode(); nodeit != mesh->GetMesh()->EndNode(); ++nodeit) 
	    {
            if(nodeit->GetStatus() != Element::Ghost)
            {
                double current_x = nodeit->RealArray(geo_coords_tag)[0]*180.0/M_PI;
                double current_y = nodeit->RealArray(geo_coords_tag)[1]*180.0/M_PI;

                if (current_x <= tmp_min_x)
                {
                    tmp_min_x = current_x;
                }

                if (current_x >= tmp_max_x)
                {
                    tmp_max_x = current_x;
                }

                if (current_y <= tmp_min_y)
                {
                    tmp_min_y = current_y;
                }

                if (current_y >= tmp_max_y)
                {
                    tmp_max_y = current_y;
                }
            }
        }

        min_geo_x = tmp_min_x;
        min_geo_y = tmp_min_y;
        max_geo_x = tmp_max_x;
        max_geo_y = tmp_max_y;
    }
    else if (grid_type == mesh::gridType::Cgrid) 
    {
        tmp_min_x = mesh->GetMesh()->BeginFace()->RealArray(geo_coords_tag)[0]*180.0/M_PI;
        tmp_max_x = mesh->GetMesh()->BeginFace()->RealArray(geo_coords_tag)[0]*180.0/M_PI;
        tmp_min_y = mesh->GetMesh()->BeginFace()->RealArray(geo_coords_tag)[1]*180.0/M_PI;
        tmp_max_y = mesh->GetMesh()->BeginFace()->RealArray(geo_coords_tag)[1]*180.0/M_PI;

        for(auto edgeit = mesh->GetMesh()->BeginFace(); edgeit != mesh->GetMesh()->EndFace(); ++edgeit) 
	    {
            if(edgeit->GetStatus() != Element::Ghost)
            {
                double current_x = edgeit->RealArray(geo_coords_tag)[0]*180.0/M_PI;
                double current_y = edgeit->RealArray(geo_coords_tag)[1]*180.0/M_PI;
                
                if (current_x <= tmp_min_x)
                {
                    tmp_min_x = current_x;
                }

                if (current_x >= tmp_max_x)
                {
                    tmp_max_x = current_x;
                }

                if (current_y <= tmp_min_y)
                {
                    tmp_min_y = current_y;
                }

                if (current_y >= tmp_max_y)
                {
                    tmp_max_y = current_y;
                }
            }
        }

        min_geo_x = tmp_min_x;
        min_geo_y = tmp_min_y;
        max_geo_x = tmp_max_x;
        max_geo_y = tmp_max_y;
    }
    else
    {
        SIMUG_ERR("Currently A and C grid is available for file interpolation");
    }

    BARRIER

    // Get coords from netcdf
    int retval;
    int fileid;

    int lon_id, lat_id;
    size_t x_size, y_size;

    if ((retval = nc_open(file_info.filename.c_str(), NC_NOWRITE, &fileid)))
    {
        NC_ERR(retval);
    }

    if((retval = nc_inq_dimid(fileid, file_info.lonname.c_str(), &lon_id)))
    {
        NC_ERR(retval);
    }

    if((retval = nc_inq_dimid(fileid, file_info.latname.c_str(), &lat_id)))
    {
        NC_ERR(retval);
    }
    
    if((retval = nc_inq_dimlen(fileid, lon_id, &x_size)))
    {
        NC_ERR(retval);
    }

    if((retval = nc_inq_dimlen(fileid, lat_id, &y_size)))
    {
        NC_ERR(retval);
    }
    
    double x_coords[x_size];
    double y_coords[y_size];

    int x_coords_id, y_coords_id;

    if ((retval = nc_inq_varid(fileid, file_info.lonname.c_str(), &x_coords_id)))
    {
        NC_ERR(retval);
    }
    
    if ((retval = nc_inq_varid(fileid, file_info.latname.c_str(), &y_coords_id)))
    {
        NC_ERR(retval);
    }

    if ((retval = nc_get_var_double(fileid, x_coords_id, x_coords)))
    {
        NC_ERR(retval);
    }
    
    if ((retval = nc_get_var_double(fileid, y_coords_id, y_coords)))
    {
        NC_ERR(retval);
    }

    double dx = x_coords[1] - x_coords[0];
    double dy = y_coords[1] - y_coords[0];

    // find x start count
    size_t x_start;
    size_t x_count;

    double max_x_begin = std::max(x_coords[0], min_geo_x);
    double min_x_end = std::min(x_coords[x_size -1], max_geo_x);

    if (max_x_begin > min_x_end)
    {
        x_count = 0;
        x_start = 0;
    }
    else if (max_x_begin == x_coords[0])
    {
        x_start = 0;
        x_count = std::min(size_t((max_geo_x - x_coords[0])/dx) + 1, x_size);
    }
    else if (max_x_begin == min_geo_x)
    {
        x_start = size_t((min_geo_x - x_coords[0])/dx);
        x_count = std::min((size_t)((max_geo_x - (x_coords[0] + dx*x_start))/dx) + 1, (x_size - x_start));
    }
    else
    {
        SIMUG_ERR("bad x_start x_count procedure");
    }

    // find y start count
    size_t y_start;
    size_t y_count;

    double max_y_begin = std::max(y_coords[0], min_geo_y);
    double min_y_end = std::min(y_coords[y_size -1], max_geo_y);

    if (max_y_begin > min_y_end)
    {
        y_count = 0;
        y_start = 0;
    }
    else if (max_y_begin == y_coords[0])
    {
        y_start = 0;
        y_count = std::min(size_t((max_geo_y - y_coords[0])/dy) + 1, y_size);
    }
    else if (max_y_begin == min_geo_y)
    {
        y_start = size_t((min_geo_y - y_coords[0])/dy);
        y_count = std::min((size_t)((max_geo_y - (y_coords[0] + dy*y_start))/dy) + 1, (y_size - y_start));
    }
    else
    {
        SIMUG_ERR("bad y_start y_count procedure");
    }

    // Find data id
    int data_variable1_id;
    int data_variable2_id;

    if ((retval = nc_inq_varid(fileid, nc_variable_u.c_str(), &data_variable1_id)))
    {
        NC_ERR(retval);
    }

    if ((retval = nc_inq_varid(fileid, nc_variable_v.c_str(), &data_variable2_id)))
    {
        NC_ERR(retval);
    }


    //Get scale factor 
    double scale_factor;
    if (file_info.scale_factor_name.size() != 0)
    {
        if ((retval = nc_get_att(fileid, data_variable1_id, file_info.scale_factor_name.c_str(), &scale_factor)))
        {
            NC_ERR(retval);
        }
    }
    else
    {
        scale_factor = 1.0;
    }
    
    //Get invalid value
    nc_type atttype;
    
    int invalid_value_int;
    short invalid_value_short;
    double invalid_value_double;

    if (file_info.invalid_value_name.size() != 0)
    {
        if ((retval = nc_inq_att(fileid, data_variable1_id, file_info.invalid_value_name.c_str(), &atttype, NULL)))
        {
            NC_ERR(retval);
        }

        if (atttype == NC_SHORT)
        {
            if ((retval = nc_get_att(fileid, data_variable1_id, file_info.invalid_value_name.c_str(), &invalid_value_short)))
            {
                NC_ERR(retval);
            }
        }
        else if ((atttype == NC_INT) or (atttype == NC_LONG))
        {
            if ((retval = nc_get_att(fileid, data_variable1_id, file_info.invalid_value_name.c_str(), &invalid_value_int)))
            {
                NC_ERR(retval);
            }
        }
        else if ((atttype == NC_FLOAT) or (atttype == NC_DOUBLE))
        {
            if ((retval = nc_get_att(fileid, data_variable1_id, file_info.invalid_value_name.c_str(), &invalid_value_double)))
            {
                NC_ERR(retval);
            }
        }
        else
        {
            SIMUG_ERR("att type is not int or double, can't read");
        }
    }
    else
    {
        invalid_value_int = std::nan("1");
        invalid_value_double = std::nan("1");
        invalid_value_short = std::nan("1");
    }

    //Get offset value
    double offset_value;
    if (file_info.offset_name.size() != 0)
    {
        if ((retval = nc_get_att(fileid, data_variable1_id, file_info.offset_name.c_str(), &offset_value)))
        {
            NC_ERR(retval);
        }
    }
    else
    {
        offset_value = 0.0;
    }

    // Assemble start and count array
    size_t* start;
    size_t* count;

    if (file_info.is_depth)
    {
        start = new size_t[4];
        count = new size_t[4];
        start[0] = netcdf_index;
        start[1] = 0;
        start[2] = y_start;
        start[3] = x_start;
        count[0] = 1;
        count[1] = 1;
        count[2] = y_count;
        count[3] = x_count;
    }
    else
    {
        start = new size_t[3];
        count = new size_t[3];
        start[0] = netcdf_index;
        start[1] = y_start;
        start[2] = x_start;
        count[0] = 1;
        count[1] = y_count;
        count[2] = x_count;
    }

    // Get data in parallel
    int vector_data_local_int1[y_count][x_count];
    int vector_data_local_int2[y_count][x_count];

    short vector_data_local_short1[y_count][x_count];
    short vector_data_local_short2[y_count][x_count];

    double vector_data_local_double1[y_count][x_count];
    double vector_data_local_double2[y_count][x_count];

    if((retval = nc_inq_var(fileid, data_variable1_id, NULL, &atttype, NULL, NULL, NULL)))
    {
        NC_ERR(retval)
    }

    
    if ((atttype == NC_SHORT))
    {
        if ((retval = nc_get_vara_short(fileid, data_variable1_id, start, count, &vector_data_local_short1[0][0])))
        {
            NC_ERR(retval)
        }
        
        if ((retval = nc_get_vara_short(fileid, data_variable2_id, start, count, &vector_data_local_short2[0][0])))
        {
            NC_ERR(retval)
        }
    }
    else if ((atttype == NC_INT) or (atttype == NC_LONG))
    {
        if ((retval = nc_get_vara_int(fileid, data_variable1_id, start, count, &vector_data_local_int1[0][0])))
        {
            NC_ERR(retval)
        }
        
        if ((retval = nc_get_vara_int(fileid, data_variable2_id, start, count, &vector_data_local_int2[0][0])))
        {
            NC_ERR(retval)
        }
    }
    else if ((atttype == NC_FLOAT) or (atttype == NC_DOUBLE))
    {
        if ((retval = nc_get_vara_double(fileid, data_variable1_id, start, count, &vector_data_local_double1[0][0])))
        {
            NC_ERR(retval)
        }

        if ((retval = nc_get_vara_double(fileid, data_variable2_id, start, count, &vector_data_local_double2[0][0])))
        {
            NC_ERR(retval)
        }
    }
    else
    {
        SIMUG_ERR("var type is not int or double, can't read");
    }
    
    // convert velocity components to double
    if (atttype == NC_SHORT)
    {
        for (size_t j = 0; j < y_count; ++j)
        {
            for(size_t i = 0; i < x_count; ++i)
            {
                double current_val_x, current_val_y;
                current_val_x = (vector_data_local_short1[j][i] == invalid_value_short) ? invalid_value_fill 
                                 : (double)(vector_data_local_short1[j][i])*scale_factor + offset_value;
                current_val_y = (vector_data_local_short2[j][i] == invalid_value_short) ? invalid_value_fill
                                 : (double)(vector_data_local_short2[j][i])*scale_factor + offset_value; 
                vector_data_local_double1[j][i] = current_val_x;
                vector_data_local_double2[j][i] = current_val_y;
            }
        }
    }           
    else if ((atttype == NC_INT) or (atttype == NC_LONG))
    {
        for (size_t j = 0; j < y_count; ++j)
        {
            for(size_t i = 0; i < x_count; ++i)
            {
                double current_val_x, current_val_y;
                current_val_x = (vector_data_local_int1[j][i] == invalid_value_int) ? invalid_value_fill 
                                 : (double)(vector_data_local_int1[j][i])*scale_factor + offset_value;
                current_val_y = (vector_data_local_int2[j][i] == invalid_value_int) ? invalid_value_fill
                                 : (double)(vector_data_local_int2[j][i])*scale_factor + offset_value;
                vector_data_local_double1[j][i] = current_val_x;
                vector_data_local_double2[j][i] = current_val_y;
            }
        }
    }
    else if ((atttype == NC_FLOAT) or (atttype == NC_DOUBLE))
    {
        for (size_t j = 0; j < y_count; ++j)
        {
            for(size_t i = 0; i < x_count; ++i)
            {
                 double current_val_x, current_val_y;
                current_val_x = (vector_data_local_double1[j][i] == invalid_value_double) ? invalid_value_fill 
                                     : (vector_data_local_double1[j][i])*scale_factor + offset_value;
                current_val_y = (vector_data_local_double2[j][i] == invalid_value_double) ? invalid_value_fill 
                                     : (vector_data_local_double2[j][i])*scale_factor + offset_value;
                vector_data_local_double1[j][i] = current_val_x;
                vector_data_local_double2[j][i] = current_val_y;
            }
        }
    }
    else
    {
        SIMUG_ERR("unknown value type");
    }


    // rotate vector to model directions
    double rotated_vector_x[y_count][x_count];
    double rotated_vector_y[y_count][x_count];

    for (size_t j = 0; j < y_count; ++j)
    {
        for (size_t i = 0; i < x_count; ++i)
        {
            double vec_x = vector_data_local_double1[j][i];
            double vec_y = vector_data_local_double2[j][i];
            std::vector<double> model_rot_vec = from_geo_2_model_vec<double>(vec_x, vec_y,
                                                                             x_coords[x_start + i],
                                                                             y_coords[y_start + j]);
            rotated_vector_x[j][i] = model_rot_vec[0];
            rotated_vector_y[j][i] = model_rot_vec[1];
       }
    }

    // get rectangle variables and interpolate
    if (grid_type == mesh::gridType::Agrid)
    {
        for(auto nodeit = mesh->GetMesh()->BeginNode(); nodeit != mesh->GetMesh()->EndNode(); ++nodeit) 
        {
            double local_x_start = x_coords[x_start];
            double local_y_start = y_coords[y_start];
            double local_x_end = x_coords[x_start + x_count-1];
            double local_y_end = y_coords[y_start + y_count-1];

            if(nodeit->GetStatus() != Element::Ghost)
            {
                double x = nodeit->RealArray(geo_coords_tag)[0]*180.0/M_PI;
                double y = nodeit->RealArray(geo_coords_tag)[1]*180.0/M_PI;

                if ((x < local_x_start) or
                    (x > local_x_end) or
                    (y < local_y_start) or
                    (y > local_y_end))
                {
                    nodeit->RealArray(vector_tag)[0] = no_extrapolation_fill;
                    nodeit->RealArray(vector_tag)[1] = no_extrapolation_fill;
                    continue;
                }
                else
                {
                    size_t x_prev_pos = (size_t)((x - local_x_start)/dx);
                    size_t y_prev_pos = (size_t)((y - local_y_start)/dy);

                    double datax_ld = rotated_vector_x[y_prev_pos][x_prev_pos];
                    double datax_lu = rotated_vector_x[y_prev_pos+1][x_prev_pos];
                    double datax_rd = rotated_vector_x[y_prev_pos][x_prev_pos+1];
                    double datax_ru = rotated_vector_x[y_prev_pos+1][x_prev_pos+1];

                    double datay_ld = rotated_vector_y[y_prev_pos][x_prev_pos];
                    double datay_lu = rotated_vector_y[y_prev_pos+1][x_prev_pos];
                    double datay_rd = rotated_vector_y[y_prev_pos][x_prev_pos+1];
                    double datay_ru = rotated_vector_y[y_prev_pos+1][x_prev_pos+1];

                    double xl = x_coords[x_start + x_prev_pos];
                    double xr = x_coords[x_start + x_prev_pos + 1];
                    double yd = y_coords[y_start + y_prev_pos];
                    double yu = y_coords[y_start + y_prev_pos + 1];


                    double curr_data_x = bilinear_interpolation(xl, xr, 
	    						                                yd, yu, 
	    						                                x , y,
	    						                                datax_ld, datax_lu,
                                                                datax_rd, datax_ru);

                    double curr_data_y = bilinear_interpolation(xl, xr, 
	    						                                yd, yu, 
	    						                                x , y,
	    						                                datay_ld, datay_lu,
                                                                datay_rd, datay_ru);

                    double abs_val = std::sqrt(curr_data_x*curr_data_x + curr_data_y*curr_data_y);

	    			nodeit->RealArray(vector_tag)[0] = 
                    (abs_val > max_abs_value) ? invalid_value_fill : curr_data_x;
                    nodeit->RealArray(vector_tag)[1] = 
                    (abs_val > max_abs_value) ? invalid_value_fill : curr_data_y;
                }
            }
        }
        mesh->GetMesh()->ExchangeData(vector_tag, NODE, 0);	
    }
    else if (grid_type == mesh::gridType::Cgrid)
    {
        for(auto edgeit = mesh->GetMesh()->BeginFace(); edgeit != mesh->GetMesh()->EndFace(); ++edgeit) 
        {
            double local_x_start = x_coords[x_start];
            double local_y_start = y_coords[y_start];
            double local_x_end = x_coords[x_start + x_count-1];
            double local_y_end = y_coords[y_start + y_count-1];

            if(edgeit->GetStatus() != Element::Ghost)
            {
                double x = edgeit->RealArray(geo_coords_tag)[0]*180.0/M_PI;
                double y = edgeit->RealArray(geo_coords_tag)[1]*180.0/M_PI;

                if ((x < local_x_start) or
                    (x > local_x_end) or
                    (y < local_y_start) or
                    (y > local_y_end))
                {
                    edgeit->RealArray(vector_tag)[0] = no_extrapolation_fill;
                    edgeit->RealArray(vector_tag)[1] = no_extrapolation_fill;
                    continue;
                }
                else
                {
                    size_t x_prev_pos = (size_t)((x - local_x_start)/dx);
                    size_t y_prev_pos = (size_t)((y - local_y_start)/dy);

                    double datax_ld = rotated_vector_x[y_prev_pos][x_prev_pos];
                    double datax_lu = rotated_vector_x[y_prev_pos+1][x_prev_pos];
                    double datax_rd = rotated_vector_x[y_prev_pos][x_prev_pos+1];
                    double datax_ru = rotated_vector_x[y_prev_pos+1][x_prev_pos+1];

                    double datay_ld = rotated_vector_y[y_prev_pos][x_prev_pos];
                    double datay_lu = rotated_vector_y[y_prev_pos+1][x_prev_pos];
                    double datay_rd = rotated_vector_y[y_prev_pos][x_prev_pos+1];
                    double datay_ru = rotated_vector_y[y_prev_pos+1][x_prev_pos+1];

                    double xl = x_coords[x_start + x_prev_pos];
                    double xr = x_coords[x_start + x_prev_pos + 1];
                    double yd = y_coords[y_start + y_prev_pos];
                    double yu = y_coords[y_start + y_prev_pos + 1];


                    double curr_data_x = bilinear_interpolation(xl, xr, 
	    						                                yd, yu, 
	    						                                x , y,
	    						                                datax_ld, datax_lu,
                                                                datax_rd, datax_ru);

                    double curr_data_y = bilinear_interpolation(xl, xr, 
	    						                                yd, yu, 
	    						                                x , y,
	    						                                datay_ld, datay_lu,
                                                                datay_rd, datay_ru);

                    double abs_val = std::sqrt(curr_data_x*curr_data_x + curr_data_y*curr_data_y);

	    			edgeit->RealArray(vector_tag)[0] = 
                    (abs_val > max_abs_value) ? invalid_value_fill : curr_data_x;
                    edgeit->RealArray(vector_tag)[1] = 
                    (abs_val > max_abs_value) ? invalid_value_fill : curr_data_y;
                }
            }
        }
        mesh->GetMesh()->ExchangeData(vector_tag, FACE, 0);
    }
    BARRIER

    // Close file
    if ((retval = nc_close(fileid)))
    {
        NC_ERR(retval);
    }
}

}