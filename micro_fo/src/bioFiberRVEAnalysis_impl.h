namespace bio
{
  template <typename I>
  void applyRVEBC(I bnd_bgn, I bnd_end, apf::Numbering * nm, las::LasOps * ops, las::Mat * k, las::Vec * f)
  {
    apf::Mesh * msh = apf::getMesh(apf::getField(nm));
    for(auto ent = bnd_bgn; ent != bnd_end; ++ent)
    {
      apf::NewArray<int> fxd_row_ids;
      int nedofs = apf::getElementNumbers(nm,*ent,fxd_row_ids);
      apf::DynamicVector zeroes(nedofs);
      zeroes.zero();
      ops->set(f,nedofs,&fxd_row_ids[0],&zeroes[0]);
      apf::Adjacent adj_vrts; // assuming only verts have nodes
      apf::getBridgeAdjacent(msh,*ent,1,0,adj_vrts);
      apf::DynamicMatrix eye(nedofs,nedofs);
      eye.zero();
      for(auto adj_vrt = adj_vrts.begin(); adj_vrt != adj_vrts.end(); ++adj_vrt)
      {
        apf::NewArray<int> fxd_col_ids;
        apf::getElementNumbers(nm,*adj_vrt,fxd_col_ids);
        ops->set(k,nedofs,&fxd_row_ids[0],nedofs,&fxd_col_ids[0],&eye(0,0));
      }
      for(int ii = 0; ii < nedofs; ++ii)
        for(int jj = 0; jj < nedofs; ++jj)
          eye(ii,jj) = ii == jj ? 1.0 : 0.0;
      ops->set(k,nedofs,&fxd_row_ids[0],nedofs,&fxd_row_ids[0],&eye(0,0));
      for(int ii = 0; ii < 3; ++ii)
        apf::fix(nm,*ent,0,ii,true);
    }
  }
  template <typename I>
  void freeRVEBC(I bnd_bgn, I bnd_end, apf::Numbering * num)
  {
    for(auto ent = bnd_bgn; ent != bnd_end; ++ent)
      for(int ii = 0; ii < 3; ++ii)
        apf::fix(num,*ent,0,ii,false);
  }
}
