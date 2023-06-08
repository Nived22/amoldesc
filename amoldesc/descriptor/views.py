from django.shortcuts import render, redirect
from . import values

# Create your views here.


def descriptor(req) :
    sm=req.GET.get('ele')
    if sm==None:
        val = values.get_properties_and_analysis('N(=O)[O]')
        return render(req, 'application/descriptor.html',context=val)
    val = values.get_properties_and_analysis(sm)
    return render(req, 'application/descriptor.html',context=val)

def search(req):
    sm=req.GET.get('ele')
    op=req.GET.get('options')
    if sm==None:
        return render(req, 'landingpages/search.html')
    if op=='fingerprint':
        return redirect('fingerprint/?ele='+sm)
    if op=='descriptor':
        return redirect('descriptor/?ele='+sm)
    return render(req, 'search.html')
