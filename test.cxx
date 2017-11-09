


	for(int i=0;i<N;i++){
		int LL = 1;
		for(int j=0;j<d;j++){
			int id=(i/LL)%L;
			int b = d*i+j;

			bst[2*b]=i;
			bst[2*b+1]=i+((id+1)%L-id)*LL;
			LL*=L;
		}
	}
