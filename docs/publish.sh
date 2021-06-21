#!/bin/bash
cd ~/.julia/dev
if [ -d "ACE_gh" ]
then
    echo "ACE_gh already exists."
else
    echo "ACE_gh does not yet exists - cloneing from git repo."
    git clone https://github.com/ACEsuit/ACEtb.jl.git ACE_gh
fi
cd ACE_gh
echo "make sure gh-pages is checked out"
git checkout gh-pages
echo "obtain doc build from ~/.julia/ACEtb/doc/build"
cp -R ../ACEtb/docs/build/* ./dev/
echo "commit changes"
git add dev/*
git commit -a -m "update online docs"
echo "push to gh-pages"
git push origin gh-pages
cd ~/.julia/dev/ACEtb/docs

