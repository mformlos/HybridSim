#include "Particle.h" 

Particle::Particle() : 
    Position { Vector3d::Zero() }, 
    Velocity { Vector3d::Zero() }, 
    Mass { 1.0 } { }
    
Particle::Particle(double aMass) : 
    Position { Vector3d::Zero() }, 
    Velocity { Vector3d::Zero() }, 
    Mass { aMass } { }
    
Particle::Particle(Vector3d aPosition, Vector3d aVelocity, double aMass) :
	Position { aPosition },
	Velocity { aVelocity },
	Mass { aMass } { }

MPCParticle::MPCParticle() : 
    Particle(), 
    CellIndex {} { }
    
MPCParticle::MPCParticle(double aMass) : 
    Particle(aMass), 
    CellIndex {} { }
    
MPCParticle::MPCParticle(Vector3d aPosition, Vector3d aVelocity, double aMass) :
    Particle(aPosition, aVelocity, aMass),
    CellIndex {} { }

    
bool MPCParticle::operator< (const MPCParticle& part) const {
    return (CellIndex < part.CellIndex); 
}
	

MDParticle::MDParticle(int Index) : 
    MPCParticle(), 
    Identifier{ Index },
    Force { Vector3d::Zero() }, 
    VerletPosition { Vector3d::Zero() }, 
    VerletList { }, 
    Bonds { } { }


MDParticle::MDParticle(int Index, double aMass) : 
    MPCParticle(aMass), 
    Identifier{ Index },
    Force { Vector3d::Zero() }, 
    VerletPosition { Vector3d::Zero() }, 
    VerletList { }, 
    Bonds { } { }

    
MDParticle::MDParticle(int Index, Vector3d aPosition, Vector3d aVelocity, double aMass) : 
    MPCParticle(aPosition, aVelocity, aMass), 
    Identifier{ Index },
    Force { Vector3d::Zero() }, 
    VerletPosition { Vector3d::Zero() }, 
    VerletList { },
    Bonds { } { }
    
    
void MDParticle::setBond(MDParticle& bonded) {
    Bonds.push_front(&bonded); 
}
    
