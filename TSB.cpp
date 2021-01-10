// HEADER

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
#include <conio.h>

using namespace std;


// Independent Functions

vector<vector<float>> initMatZeros( int a, int b ) {
	// Initializes an 'a' by 'b' matrix
	vector<vector<float>> matrix;
	vector<float> line;
	float zeroval = 0;
	
	// Get line, each line is the same
	for (int i = 0; i < b; ++i) {
		line.push_back(zeroval);
	}
	
	// Concatenate lines up to a
	for (int i = 0; i < a; ++i) {
		matrix.push_back(line);
	}
	
	return matrix;
}


vector<float> concvect(vector<float> v1, vector<float> v2 ) {
	// Concatenate vector one with vector 2
	vector<float> outv = v1;
	
	for (int i = 0; i < v2.size(); ++i) {
		outv.push_back( v2[i] );
	}
	
	return outv;
}


vector<float> vectScal(vector<float> v, float s ) {
	// Multiply vector with scalar
	
	for (int i = 0; i < v.size(); ++i) {
		v[i] = v[i] * s;
	}
	
	return v;
}


vector<vector<float>> matScal( vector<vector<float>> m, float s ) {
	// Product of matrix with scalar
	vector<vector<float>> mout = initMatZeros( m.size(), m[0].size() );
	
	for (int i = 0; i < m.size(); ++i) {
		for (int j = 0; j < m[0].size(); ++j) {
			mout[i][j] = m[i][j] * s;
		}
	}
	
	return mout;
}


vector<vector<float>> matTr( vector<vector<float>> m ) {
	// Transpose matrix
	vector<vector<float>> outm = initMatZeros( m[0].size(), m.size() );
	
	// Add first elements in matrix
	for (int i = 0; i < outm.size(); ++i) {
		for (int j = 0; j < outm[0].size(); ++j) {
			outm[i][j] = m[j][i];
		}
	}
	
	return m;
}


vector<vector<float>> matProd( vector<vector<float>> m1, vector<vector<float>> m2 ) {
	// Matrix multiplication of m1 and m2
	vector<vector<float>> outm = initMatZeros( m1.size(), m2[0].size() );
	
	for (int i = 0; i < outm.size(); ++i ) {
		for (int j = 0; j < outm[0].size(); ++j ) {
			for (int k = 0; k < m2.size(); ++k ) { 
				outm[i][j] += m1[i][k] * m2[k][j];
			}
		}
	}
	
	return outm;
}


vector<vector<float>> outProd( vector<float> v0, vector<float> v1) {
	// Outer product of two vectors
	vector<vector<float>> m = initMatZeros( v0.size(), v0.size() );
	
	for (int i = 0; i < v0.size(); ++i) {
		for (int j = 0; j < v0.size(); ++j) {
			m[i][j] = v0[i] * v1[j];
		}
	}
	
	return m;
}


vector<float> subvector(vector<float> invect, int a, int b) {
	// Return a subvector of invect, from indices a to b
	vector<float> outvect;
	float data;
	
	for (int i = a; i <= b; ++i) {
		data = invect[i];
		outvect.push_back(data);
	}
	
	return outvect;
}


vector<vector<float>> cutmatrix(vector<vector<float>> m, int c0, int c1) {
	// Cut lines c0 to c1 and collumns c0 to c1 from square matrix m
	vector<vector<float>> mout = m;
	
	// Cut lines
	mout.erase( mout.begin() + c0, mout.begin() + c1);
	
	// Cut collumn c
	for (int i = 0; i < mout.size(); ++i) {
		mout[i].erase( mout[i].begin() + c0, mout[i].begin() + c1);
	}
	
	return mout;
}


vector<float> addvector(vector<float> v0, vector<float> v1 ) {
	// Adds v0 and v1 element-wise
	vector<float> vout = v0;
	
	for (int i = 0; i < v0.size(); ++i) {
		vout[i] = v0[i] + v1[i];
	}
	
	return vout;
}


vector<vector<float>> addMat( vector<vector<float>> m0, vector<vector<float>> m1 ) {
	// Add m0 and m1 element-wise
	vector<vector<float>> mout = initMatZeros( m0.size(), m0[0].size() );
	
	for (int i = 0; i < m0.size(); ++i) {
		for (int j = 0; j < m0[0].size(); ++j) {
			mout[i][j] = m0[i][j] + m1[i][j];
		}
	}
	
	return mout;
}


vector<float> solvelin( vector<vector<float>> A, vector<float> b ) {
	// Solve linear system A*x = b via gauss-jordan method
	vector<float> x = b; // init x as b
	bool isinvertible = true;
	float piv;
	
	// Loop over all collumns to produce an upper triangular system
	for (int j = 0; j < A[0].size(); ++j) {
		// Check if A[j][j] is 0. if it is, permute line i with next line that isn't 0
		if (A[j][j] == 0) {
			vector<float> linej = A[j];
			float valj = x[j];
			bool foundnonzero = false;
			for (int i = j; i < A.size(); ++i) {
				if (A[i][j] != 0) {
					A[j] = A[i];
					A[i] = linej;
					x[j] = x[i];
					x[i] = valj;
					i+= A.size();
					foundnonzero = true;
				}
			}
			
			if (!foundnonzero) {
				isinvertible = false;
			}
		}
		
		// Loop over all the lines greater than j
		for (int i = j+1; i < A.size(); ++i) {
			// Check if A[i][j] is 0. if it isn't, eliminate
			if (A[i][j] != 0) {
				piv = -1 * A[i][j] / A[j][j];
				A[i] = addvector( A[i], vectScal(A[j], piv) );
				x[i] += x[j] * piv;
			}
		}
	}
	
	if (isinvertible) {
		// Redo the thing but like in reverse
		for (int j = A[0].size()-1; j >= 0; --j) {
			// Assume that previous step was a sucess. No need to check pivots == 0
			// Loop over all the lines smaller than j
			for (int i = j-1; i >= 0; --i) {
				// Check if A[i][j] is 0. if it isn't, eliminate
				if (A[i][j] != 0) {
					piv = -1 * A[i][j] / A[j][j];
					A[i] = addvector( A[i], vectScal(A[j], piv) );
					x[i] += x[j] * piv;
				}
			}
		}
		
		// Divide each value of x with its corresponding pivot
		for (int i = 0; i < x.size(); ++i) {
			x[i] = x[i] / A[i][i];
		}
	} else {
		x = {0, 1, 1, 0, 0, 1, 1, 0};
	}
	
	return x;
}


vector<float> readline(string stringline) {
	// Read a string line of values separated by spaces
	vector<float> linedat;
	int lasti = 0;
	string thischar;
	string thisval;
	const char * c;
	
	// Read up to last value
	for (int i = 0; i < stringline.length(); ++i) {
		// Get character i of string
		thischar = stringline.at(i);
		
		// Check if character i is a space
		if (thischar.compare(" ") == 0) {
			// If it's a space, from lasti to i is a value. push back value to linedat
			thisval = stringline.substr(lasti, i);
			if (thisval.compare(" ; ") != 0) {
				c = thisval.c_str();
				linedat.push_back(std::atof(c));
			}
			// Change lasti to be i
			lasti = i;
		}
	}
	
	// Read last value into linedat
	thisval = stringline.substr(lasti, stringline.length()-1);
	c = thisval.c_str();
	linedat.push_back(std::atof(c));
	
	return linedat;
}


// Classes

class frameClass {
	public:
		// Properties
		float E;
		float d;
		vector<vector<float>> v;
		vector<vector<float>> f;
		vector<float> fc;
		vector<int> f_adress;
		vector<int> c;
		vector<float> e;
		vector<vector<float>> e_adress;
		
		// Methods
		vector<float> calcStress() {
			// Get unit components of each edges
			vector<vector<float>> unitedges;
			vector<float> thisedge;
			vector<float> edgelen;
			float mag;
			float thisval;
			int v0, v1, h, t;
			for (int i = 0; i < this->e.size(); ++i) {
				// Get start and end vertices of this edge
				v0 = (int)this->e_adress[i][0];
				v1 = (int)this->e_adress[i][1];
				thisedge = {};
				for (int j = 0; j < this->v[v0].size(); ++j) {
					thisedge.push_back( this->v[v1][j] - this->v[v0][j] );
				}
				// Normalize vector
				mag = 0;
				for (int j = 0; j < thisedge.size(); ++j) {
					mag += thisedge[j] * thisedge[j];
				}
				mag = sqrt(mag);
				edgelen.push_back(mag); // Save length of edge i
				mag = 1/mag;
				for (int j = 0; j < thisedge.size(); ++j) {
					thisedge[j] = thisedge[j] * mag;
				}
				// Add edge to edge list
				unitedges.push_back( thisedge );
			}
			
			// Initialize stiffness matrix
			vector<vector<float>> stiffmat = initMatZeros( 3*this->v.size(), 3*this->v.size() );
			
			// Loop over all the edges to add to the stiffness matrix
			for (int n = 0; n < unitedges.size(); ++n) {
				// Float vector kn0 to build kn
				vector<float> kn0 = concvect( vectScal(unitedges[n], -1) , unitedges[n] );
				
				// Matrix multiply kn0 with itself transposed
				vector<vector<float>> kn = outProd( kn0, kn0 );
				// Multiply kn with scalar 
				kn = matScal( kn, this->E * this->e[n] / edgelen[n] );
				
				// Get start and end vertices of this edge
				v0 = (int)this->e_adress[n][0];
				v1 = (int)this->e_adress[n][1];
				
				for (int i = 0; i < kn.size(); ++i) {
					for (int j = 0; j < kn[0].size(); ++j) {
						if (i < 3) {
							h = 3*(v0+1) - 3 + i;
						} else {
							h = 3*(v1+1) - 6 + i;
						}
						if (j < 3) {
							t = 3*(v0+1) - 3 + j;
						} else {
							t = 3*(v1+1) - 6 + j;
						}
						stiffmat[h][t] += kn[i][j];
					}
				}
			}
			
			// Add edges weights to each vertex z forces
			vector<float> fc2 = this->fc;
			float this_weight;
			for (int i = 0; i < this->e.size(); ++i) {
				// Get start and end vertices of this edge
				v0 = (int)this->e_adress[i][0];
				v1 = (int)this->e_adress[i][1];
				
				// Get this edge weight, halved
				this_weight = this->e[i] * edgelen[i] * this->d * 4.905;
				
				// Add the half weight to v0 and v1 forces
				fc2[3*v0 + 2] -= this_weight;
				fc2[3*v1 + 2] -= this_weight;
			}
			
			// Remove constrained nodes from system matrix and forces vector
			vector<int> c2 = this->c;
			for (int i = 0; i < c2.size(); ++i) {
				stiffmat = cutmatrix(stiffmat, 3*c2[i], 3*c2[i]+3 );
				// Lower index of vertex bigger than the one removed
				for (int j = i; j < c2.size(); ++j) {
					if (c2[i] < c2[j]) {
						c2[j] -= 1;
					}
				}
				fc2.erase( fc2.begin() + 3*c2[i], fc2.begin() + 3*c2[i] + 3 );
			}
			
			// Solve linear system stiffmat * u = fc2 for edge displacement
			vector<float> usol = solvelin(stiffmat, fc2);
			
			// Check if usol means success of fail of linear solve
			bool isinv;
			if (usol.size() == 8) {
				isinv = false;
				vector<float> check_inv = {0, 1, 1, 0, 0, 1, 1, 0};
				for (int i = 0; i < usol.size(); ++i) {
					if (usol[i] != check_inv[i]) {
						isinv = true;
					}
				}
			} else {
				isinv = true;
			}
			
			// Declare stress vector output and check if isinv
			vector<float> stress;
			if (isinv) {
				// Convert usol into matrix U with all vertex
				vector<vector<float>> U;
				vector<float> zerv = { 0, 0, 0 };
				int k = 0;
				for (int i = 0; i < this->v.size(); ++i) {
					// Check if vertex i is one of the constrained ones
					int i_in_c = 0;
					for (int j = 0; j < this->c.size(); ++j) {
						if (i == this->c[j]) {
							if (i_in_c == 0) {
								i_in_c = 1;
								U.push_back(zerv);
								k = k + 1;
							}
						}
					}
					if (i_in_c == 0) {
						vector<float> thisu = {};
						for (int j = 0; j < 3; ++j) {
							thisu.push_back( usol[3*(i-k)+j] );
						}
						U.push_back(thisu);
					}
				}
				
				// Compute new vertex position
				vector<vector<float>> new_v = addMat( U, this->v );
				
				// Use U to compute new edges length
				vector<float> newedgelen;
				for (int i = 0; i < this->e.size(); ++i) {
					// Get start and end vertices of this edge
					v0 = (int)this->e_adress[i][0];
					v1 = (int)this->e_adress[i][1];
					thisedge = {};
					for (int j = 0; j < new_v[v0].size(); ++j) {
						thisedge.push_back( new_v[v1][j] - new_v[v0][j] );
					}
					// Normalize vector
					mag = 0;
					for (int j = 0; j < thisedge.size(); ++j) {
						mag += thisedge[j] * thisedge[j];
					}
					mag = sqrt(mag);
					newedgelen.push_back(mag); // Save length of edge i
				}
				
				// Compute vector elongation and stress
				vector<float> elong = edgelen;
				stress = edgelen;
				for (int i = 0; i < elong.size(); ++i) {
					elong[i] = ( newedgelen[i] - edgelen[i] ) / edgelen[i];
					stress[i] = elong[i] * this->E;
				}
			} else {
				stress = {0, 1, 1, 0, 0, 1, 1, 0};
			}
			
			// Finally, output
			return stress;
		}
};


// Frame creation function

frameClass readFrame(string inpfilename) {
	// Create frame object
	frameClass frame;
	
	// Read frame input file
	std::ifstream file(inpfilename);
    string str; 
	string firstchar;
	string linedata;
	const char * c;
	vector<int> numreads = {0, 0, 0, 0, 0, 0};
	vector<string> props = {"E", "d", "c", "v", "f", "e"};
    while (std::getline(file, str)) {
        // Process each str line
		
		// Check first character of string
		if (str.length() > 0) {
			firstchar = str.at(0);
			linedata = str.substr(3);
			if (firstchar.compare("E") == 0) {
				c = linedata.c_str();
				frame.E = std::atof(c);
				numreads[0] += 1;
			} else if (firstchar.compare("d") == 0) {
				c = linedata.c_str();
				frame.d = std::atof(c);
				numreads[1] += 1;
			} else if (firstchar.compare("c") == 0) {
				c = linedata.c_str();
				frame.c.push_back(std::atof(c));
				numreads[2] += 1;
			} else if (firstchar.compare("v") == 0) {
				vector<float> vectdata = readline(linedata);
				frame.v.push_back(vectdata);
				numreads[3] += 1;
			} else if (firstchar.compare("f") == 0) {
				vector<float> vectdata = readline(linedata);
				frame.f_adress.push_back(vectdata[0]);
				frame.f.push_back( subvector(vectdata, 1, vectdata.size()-1) );
				numreads[4] += 1;
			} else if (firstchar.compare("e") == 0) {
				vector<float> vectdata = readline(linedata);
				frame.e.push_back(vectdata[0]);
				frame.e_adress.push_back( subvector(vectdata, 2, vectdata.size()-1) );
				numreads[5] += 1;
			}
		}
    }
	
	// Check for empty parameters
	for (int i = 0; i < numreads.size(); ++i) {
		if (numreads[i] == 0) {
			cout << "Warning: property " << props[i] << " has not been found in frame file." << endl;
			cout << endl;
		}
	}
	
	// Create fc vector
	vector<float> f2add;
	for (int i = 0; i < frame.v.size(); ++i) {
		f2add = {0, 0, 0};
		for (int j = 0; j < frame.f_adress.size(); ++j) {
			if (i == frame.f_adress[j]) {
				f2add = addvector( f2add, frame.f[j] );
			}
		}
		
		for (int j = 0; j < f2add.size(); ++j) {
			frame.fc.push_back( f2add[j] );
		}
	}
	
	return frame;
}


// Write data to file function

int writeVectText( string outfilename, vector<float> data ) {
	// Write data in data to csv format
	
	ofstream outfile;
	outfile.open( outfilename );
	for (int i = 0; i < data.size(); ++i) {
		outfile << data[i] << " \n";
	}
	outfile.close();
	
	return 0;
}


// Main

int main(int argc, char* argv[]) {
	// Define Strings
	string inpfilename;
	string outfilename;
	
	// Print name of program
	cout << "" << endl;
	cout << "Truss Stress [version 0.1 - private release]" << endl;
	cout << "(c) 2021 Alexis Angers [https://github.com/Alexsimulation]. Private and educational use only." << endl;
	cout << "" << endl;
	
	// Input management
	if (argc > 1) {
		inpfilename = argv[1];
		cout << " Input file : " << inpfilename << endl;
		if (argc > 2) {
			outfilename = argv[2];
		} else {
			outfilename = "none";
		}
		cout << "Output file : " << outfilename << endl;
		cout << endl;
		
		// Read frame file
		frameClass frame = readFrame(inpfilename);
		
		// Compute stress from frame
		vector<float> stress = frame.calcStress();
		
		// Check if stress means success of fail of linear solve
		bool isinv;
		if (stress.size() == 8) {
			isinv = false;
			vector<float> check_inv = {0, 1, 1, 0, 0, 1, 1, 0};
			for (int i = 0; i < stress.size(); ++i) {
				if (stress[i] != check_inv[i]) {
					isinv = true;
				}
			}
		} else {
			isinv = true;
		}
		
		if (isinv) {
			// Output stress to command wondow
			for (int i = 0; i < stress.size(); ++i) {
				cout << "e" << i << " : " << stress[i] << " Pa" << endl;
			}
			
			// Write stress to output file
			if (outfilename.compare("none") != 0) {
				writeVectText( outfilename, stress );
				cout << endl << "Data written to output file" << endl;
			}
		} else {
			cout << "Error: non-invertible stiffness matrix. Verify that the truss has at least 3 constrained vertices." << endl;
		}
	} else {
		cout << "No input file. Program will exit now." << endl;
	}
	
	return 0;
}



