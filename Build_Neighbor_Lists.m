function [FirstNeighbors,SecondNeighbors,DistantAtoms,AtomData,iA,jA] = Build_Neighbor_Lists(AL)

%This function runs over all atoms and for each defines a list of first
%neighbor, second neighbor and distant atoms. Also, based on the output of
%Build_Basis, it will assign the nb limits to each atom (the first and last basis
%functions involving that atom, and the mu limits (the first and last CGTF
%indexes involving that atom).
%The difference between nb and mu is that nb = nb + 1 for s,p,d,f, shells, but
%mu = mu + size_shell (currently 1,3,6,10, but sometime to be contracted to
%1,3,5,7,etc)

NAtom = size(AL,1);
t = 1;
s = 1;
r = 1;
%The thresholds for defining first, second and distant neighbors should be
%adjusted based on numerical experimentation. Also it would be possible to
%define a threshold list so that some atoms can be treated more accurately,
%or the threshold could vary for heavier or lighter atoms.

thr1 = 5; %bohr
thr2 = 8; %bohr

FirstNeighbors = cell(NAtom,1);
SecondNeighbors = cell(NAtom,1);
DistantAtoms = cell(NAtom,1);
AtomData = [];

for iA = 1:NAtom
	for jA = iA+1:NAtom
		if (norm(AL(jA,:)-AL(iA,:)) <= thr1) %example 3 A
			FirstNeighbors{iA}(t) = jA;
			t = t+1;
		elseif (norm(AL(jA,:)-AL(iA,:)) > thr1 && norm(AL(jA,:)-AL(iA,:)) <= thr2) %example 6 A
			SecondNeighbors{iA}(s) = jA;
			s = s + 1;
		else
			DistantAtoms{iA}(r) = jA;
			r = r + 1;
		end
	end
end


end