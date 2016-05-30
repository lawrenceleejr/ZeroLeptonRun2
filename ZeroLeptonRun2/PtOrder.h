#ifndef PtOrder_h
#define PtOrder_h

template <class Proxy>
static bool PtOrder(const Proxy& first, const Proxy& second) { return first.Pt() > second.Pt();}

#endif
