%{
---------------------------------------------------------------------------
Author: Yu-Huan Wang 
    (Kim Lab at UIUC) - yuhuanw2@illinois.edu
    Creation date: 3/19/2024
    Last updated at 3/19/2024

Description: this script sort the files by natural order, it only works for
fileName containing 't-time' (e.g. t30, t60) 

---------------------------------------------------------------------------
%}    

function newList = sortNat( list)

    timeP = string( regexp( {list.name}, 't\d+', 'match', 'once'))'; % e.g. Cy3_spotsMesh_t30.mat
    [~, timeOrder] = sort( str2double( erase( timeP, "t")));
    newList = list( timeOrder);
end