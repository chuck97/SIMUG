#include "fct_zalesak.hpp"

using namespace std;
using namespace INMOST;

namespace SIMUG
{

void FCT_Procedure(INMOST::Mesh* mesh,
                   const std::vector<std::vector<std::vector<double>>>& M_C_minus_M_L,
                   INMOST::Tag scal_tag,
                   INMOST::Tag scal_low_tag,
                   INMOST::Tag scal_high_tag,
                   INMOST::Tag node_id_tag,
                   double fct_cd)
{
    // calculate fluxes
    std::vector<std::vector<double>> Fe_vector(M_C_minus_M_L.size());
    int tr_num = 0;

    for(auto trianit = mesh->BeginCell();
             trianit != mesh->EndCell();
             ++trianit)
    {
        if (trianit->GetStatus() != Element::Ghost)
        {
            INMOST::ElementArray<Node> local_nodes = trianit->getNodes();
            
            std::vector<double> scal_n_vec = 
            {
                local_nodes[0]->Real(scal_tag),
                local_nodes[1]->Real(scal_tag),
                local_nodes[2]->Real(scal_tag)
            };

            std::vector<double> scal_H_vec = 
            {
                local_nodes[0]->Real(scal_high_tag),
                local_nodes[1]->Real(scal_high_tag),
                local_nodes[2]->Real(scal_high_tag)
            };

            std::vector<double> vec = (fct_cd - 1.0)*scal_n_vec + scal_H_vec;

            Fe_vector[tr_num] = M_C_minus_M_L[tr_num]*vec;

            ++tr_num;
        }
    }
    BARRIER

    // make fluxes tags on triangles
    INMOST::Tag Fe_trian_tag;
    tr_num = 0;
    Fe_trian_tag = mesh->CreateTag("Fe_trian", DATA_REAL, CELL, NONE, 3);
    mesh->SetFileOption("Tag:Fe_trian", "nosave");
    
    for(auto trianit = mesh->BeginCell();
             trianit != mesh->EndCell();
             ++trianit)
    {
        if (trianit->GetStatus() != Element::Ghost)
        {
            trianit->RealArray(Fe_trian_tag)[0] = Fe_vector[tr_num][0];
            trianit->RealArray(Fe_trian_tag)[1] = Fe_vector[tr_num][1];
            trianit->RealArray(Fe_trian_tag)[2] = Fe_vector[tr_num][2];
            ++tr_num;
        }
    }

    // exchange data
    mesh->ExchangeData(Fe_trian_tag, CELL, 0);
    BARRIER

    // make alpha tags on triangles
    INMOST::Tag alpha_trian_tag;
    tr_num = 0;
    alpha_trian_tag = mesh->CreateTag("alpha_trian", DATA_REAL, CELL, NONE, 1);

    // calculate coefficients
    CalculateAlpha(mesh, Fe_trian_tag, alpha_trian_tag, scal_tag, scal_low_tag, node_id_tag);

    // initialize scal with scal_low
    for(auto nodeit = mesh->BeginNode();
             nodeit != mesh->EndNode();
             ++nodeit)
    {
        if(nodeit->GetStatus() != Element::Ghost)
        {
            nodeit->Real(scal_tag) = nodeit->Real(scal_low_tag);
        }
    }
    
    // exchange data
    mesh->ExchangeData(scal_tag, NODE, 0);
    BARRIER

    // update scalar value
    for(auto trianit = mesh->BeginCell();
             trianit != mesh->EndCell();
             ++trianit)
    {
        INMOST::ElementArray<Node> local_nodes = trianit->getNodes();

        for (int i = 0; i < 3; ++i)
        {
            if (local_nodes[i]->GetStatus() != Element::Ghost)
            {
                local_nodes[i].Real(scal_tag) += trianit->RealArray(Fe_trian_tag)[i]*
                                                 trianit->Real(alpha_trian_tag);
            }
        }
    }

    mesh->ExchangeData(scal_tag, NODE, 0);
    BARRIER

    // Delete temporal Tags
    mesh->DeleteTag(Fe_trian_tag, CELL);
    mesh->DeleteTag(alpha_trian_tag, CELL);
}

void CalculateAlpha(INMOST::Mesh* mesh,
                    INMOST::Tag Fe_trian_tag,
                    INMOST::Tag alpha_trian_tag,
                    INMOST::Tag scal_tag,
                    INMOST::Tag scal_low_tag,
                    INMOST::Tag node_id_tag)
{
    // calculate scal_node_str: 0 - min, 1 - max
    INMOST::Tag scal_node_str_tag;
    scal_node_str_tag = mesh->CreateTag("scal_node_str", DATA_REAL, NODE, NONE, 2);
    mesh->SetFileOption("Tag:scal_node_str", "nosave");
    
    for(auto nodeit = mesh->BeginNode();
             nodeit != mesh->EndNode();
             ++nodeit)
    {
        if(nodeit->GetStatus() != Element::Ghost)
        {
            double scal_low = nodeit->Real(scal_low_tag);
            double scal_n = nodeit->Real(scal_tag);
            nodeit->RealArray(scal_node_str_tag)[0] = std::min(scal_low, scal_n);
            nodeit->RealArray(scal_node_str_tag)[1] = std::max(scal_low, scal_n);
        }
    }
    mesh->ExchangeData(scal_node_str_tag, NODE, 0);
    BARRIER

    // calculate scal_trian_str_str: 0 - min, 1 - max
    INMOST::Tag scal_trian_str_str_tag;
    scal_trian_str_str_tag = mesh->CreateTag("scal_trian_str_str_tag", DATA_REAL, CELL, NONE, 2);
    mesh->SetFileOption("Tag:scal_trian_str_str_tag", "nosave");
    
    for(auto trianit = mesh->BeginCell();
             trianit != mesh->EndCell();
             ++trianit)
    {
        if(trianit->GetStatus() != Element::Ghost)
        {
            INMOST::ElementArray<Node> local_nodes = trianit->getNodes();
            double min1 = local_nodes[0]->RealArray(scal_node_str_tag)[0];
            double min2 = local_nodes[1]->RealArray(scal_node_str_tag)[0];
            double min3 = local_nodes[2]->RealArray(scal_node_str_tag)[0];

            double max1 = local_nodes[0]->RealArray(scal_node_str_tag)[1];
            double max2 = local_nodes[1]->RealArray(scal_node_str_tag)[1];
            double max3 = local_nodes[2]->RealArray(scal_node_str_tag)[1];

            double min_trian = std::min(std::min(min1, min2), min3);
            double max_trian = std::max(std::max(max1, max2), max3);

            trianit->RealArray(scal_trian_str_str_tag)[0] = min_trian;
            trianit->RealArray(scal_trian_str_str_tag)[1] = max_trian;
        }
    }
    mesh->ExchangeData(scal_trian_str_str_tag, CELL, 0);
    BARRIER

    // calculate m_node_min_max: 0 - min, 1 - max
    INMOST::Tag scal_node_min_max_tag;
    scal_node_min_max_tag = mesh->CreateTag("m_node_min_max", DATA_REAL, NODE, NONE, 2);
    mesh->SetFileOption("Tag:m_node_min_max", "nosave");
    
    for(auto nodeit = mesh->BeginNode();
             nodeit != mesh->EndNode();
             ++nodeit)
    {
        if(nodeit->GetStatus() != Element::Ghost)
        {
            INMOST::ElementArray<Cell> local_trians = nodeit->getCells();
            double node_min = std::numeric_limits<double>::max();;
            double node_max = std::numeric_limits<double>::min();;
            for (size_t i = 0; i < local_trians.size(); ++i)
            {
                if (local_trians[i]->RealArray(scal_trian_str_str_tag)[0] < node_min)
                {
                    node_min = local_trians[i]->RealArray(scal_trian_str_str_tag)[0];
                }

                if (local_trians[i]->RealArray(scal_trian_str_str_tag)[1] > node_max)
                {
                    node_max = local_trians[i]->RealArray(scal_trian_str_str_tag)[1];
                }
            }
            nodeit->RealArray(scal_node_min_max_tag)[0] = node_min;
            nodeit->RealArray(scal_node_min_max_tag)[1] = node_max;
        }
    }
    mesh->ExchangeData(scal_node_min_max_tag, NODE, 0);
    BARRIER

    // calculate q_node_min_max: 0 - min, 1 - max
    INMOST::Tag q_node_min_max_tag;
    q_node_min_max_tag = mesh->CreateTag("q_node_min_max", DATA_REAL, NODE, NONE, 2);
    mesh->SetFileOption("Tag:q_node_min_max_tag", "nosave");
    
    for(auto nodeit = mesh->BeginNode();
             nodeit != mesh->EndNode();
             ++nodeit)
    {
        if(nodeit->GetStatus() != Element::Ghost)
        {
            double node_scal_min = nodeit->RealArray(scal_node_min_max_tag)[0];
            double node_scal_max = nodeit->RealArray(scal_node_min_max_tag)[1];

            double scal_low = nodeit->Real(scal_low_tag);

            nodeit->RealArray(q_node_min_max_tag)[0] = node_scal_min - scal_low;
            nodeit->RealArray(q_node_min_max_tag)[1] = node_scal_max - scal_low;
        }
    }
    mesh->ExchangeData(q_node_min_max_tag, NODE, 0);
    BARRIER

    // calculate p_node_neg_pos: 0 - negative, 1 - positive
    INMOST::Tag p_node_neg_pos_tag;
    p_node_neg_pos_tag = mesh->CreateTag("p_node_neg_pos", DATA_REAL, NODE, NONE, 2);
    mesh->SetFileOption("Tag:p_node_neg_pos", "nosave");
    
    for(auto nodeit = mesh->BeginNode();
             nodeit != mesh->EndNode();
             ++nodeit)
    {
        if(nodeit->GetStatus() != Element::Ghost)
        {
            double sum_neg = 0.0;
            double sum_pos = 0.0;
            INMOST::ElementArray<Cell> local_trians = nodeit->getCells();
            for (size_t i = 0; i < local_trians.size(); ++i)
            {
                INMOST::ElementArray<Node> local_nodes = local_trians[i]->getNodes();
                for (int j = 0; j < 3; ++j)
                {
                    if (local_nodes[j]->Integer(node_id_tag) == nodeit->Integer(node_id_tag))
                    {
                        sum_neg += std::min(local_trians[i]->RealArray(Fe_trian_tag)[j], 0.0);
                        sum_pos += std::max(local_trians[i]->RealArray(Fe_trian_tag)[j], 0.0);
                    }
                }
            }

            nodeit->RealArray(p_node_neg_pos_tag)[0] = sum_neg;
            nodeit->RealArray(p_node_neg_pos_tag)[1] = sum_pos;
        }
    }
    mesh->ExchangeData(p_node_neg_pos_tag, NODE, 0);
    BARRIER

    // calculate R_node_neg_pos: 0 - negative, 1 - positive
    INMOST::Tag R_node_neg_pos_tag;
    R_node_neg_pos_tag = mesh->CreateTag("R_node_neg_pos", DATA_REAL, NODE, NONE, 2);
    mesh->SetFileOption("Tag:R_node_neg_pos", "nosave");
    
    for(auto nodeit = mesh->BeginNode();
             nodeit != mesh->EndNode();
             ++nodeit)
    {
        if(nodeit->GetStatus() != Element::Ghost)
        {
            double P_neg = nodeit->RealArray(p_node_neg_pos_tag)[0];
            double P_pos = nodeit->RealArray(p_node_neg_pos_tag)[1];

            double Q_neg = nodeit->RealArray(q_node_min_max_tag)[0];
            double Q_pos = nodeit->RealArray(q_node_min_max_tag)[1];

            if (std::abs(P_neg) < REAL_MIN_ABS_VAL)
            {
                nodeit->RealArray(R_node_neg_pos_tag)[0] = 0.0;
            }
            else
            {
                nodeit->RealArray(R_node_neg_pos_tag)[0] = std::min(1.0, Q_neg/P_neg);
            }

            if (fabs(P_pos) < REAL_MIN_ABS_VAL)
            {
                nodeit->RealArray(R_node_neg_pos_tag)[1] = 0.0;
            }
            else
            {
                nodeit->RealArray(R_node_neg_pos_tag)[1] = std::min(1.0, Q_pos/P_pos);
            }   
        }
    }

    mesh->ExchangeData(R_node_neg_pos_tag, NODE, 0);
    BARRIER

    // calculate trian_alpha
    for(auto trianit = mesh->BeginCell();
             trianit != mesh->EndCell();
             ++trianit)
    {
        if(trianit->GetStatus() != Element::Ghost)
        {
            INMOST::ElementArray<Node> local_nodes = trianit->getNodes();
            std::vector<double> Fe_vec = 
            {
                trianit->RealArray(Fe_trian_tag)[0],
                trianit->RealArray(Fe_trian_tag)[1],
                trianit->RealArray(Fe_trian_tag)[2]
            };
            
            std::vector<double> R_vec = 
            {
                (Fe_vec[0] < 0.0) ? local_nodes[0]->RealArray(R_node_neg_pos_tag)[0]: local_nodes[0]->RealArray(R_node_neg_pos_tag)[1],
                (Fe_vec[1] < 0.0) ? local_nodes[1]->RealArray(R_node_neg_pos_tag)[0]: local_nodes[1]->RealArray(R_node_neg_pos_tag)[1],
                (Fe_vec[2] < 0.0) ? local_nodes[2]->RealArray(R_node_neg_pos_tag)[0]: local_nodes[2]->RealArray(R_node_neg_pos_tag)[1]
            };

            trianit->Real(alpha_trian_tag) = std::min(std::min(R_vec[0], R_vec[1]), R_vec[2]); 
        }
    }
    mesh->ExchangeData(alpha_trian_tag, CELL, 0);
    BARRIER

    // delete all temporal tags
    mesh->DeleteTag(scal_node_str_tag, NODE);
    mesh->DeleteTag(scal_trian_str_str_tag, CELL);
    mesh->DeleteTag(scal_node_min_max_tag, NODE);
    mesh->DeleteTag(q_node_min_max_tag, NODE);
    mesh->DeleteTag(p_node_neg_pos_tag, NODE);
    mesh->DeleteTag(R_node_neg_pos_tag, NODE);
    BARRIER
}

}