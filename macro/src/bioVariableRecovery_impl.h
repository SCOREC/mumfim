namespace bio
{
  template <typename T>
    void getEntsFieldComponent(apf::Field * fld,
                               int crd,
                               std::vector<double>& crds,
                               T begin,
                               T end)
  {
    for(T it = begin; it != end; it++)
      getEntFieldComponent(fld,crd,crds,*it);
  }
}
