from django.shortcuts import render, redirect
from . import code

# Create your views here.

def fingerprint(req):
    smile= req.GET.get('ele')
    fp = code.fingerprint_generation(smile)
    return render(req,'application/fingerprint.html',context=fp)

def index(req):
    return render(req, 'landingpages/coverpage.html')

def search(req):
    sm=req.GET.get('ele')
    op=req.GET.get('options')
    if sm==None:
        return render(req, 'landingpages/search.html')
    if op=='fingerprint':
        return redirect('fingerprint/?ele='+sm)
    if op=='descriptor':
        return redirect('descriptor/?ele='+sm)
    return render(req, 'landingpages/search.html')

def contactus(req):
    return render(req, 'landingpages/contactus.html')

def developers(req):
    return render(req,'landingpages/developers.html')