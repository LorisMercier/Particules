///INFORMATIONS PROJET :
///Aimant la physique mécanique, mon projet de LIFAMI a pour but de reproduire le plus fidèlement les réactions entre des balles
///et d'autres éléments (Des parois, d'autres balles, (et bientôt de l'eau si possible)). J'essaye donc, à mon échelle, de créer
///un moteur graphique autour d'objet circulaire.


///NOUVEAUTE DE LA SEMAINE :
///Voir le README

#include <Grapic.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

using namespace grapic;
using namespace std;

const int DIMW = 850;
const int DIMH = 700;
const int MAX_PART = 40;
const int MAX_LIGNE = 11;
const int MAXR = 150;
const int SURFACE_EAU = 160;
const int FACTEUR_MASSE = 60000;


const float deltaT = 0.04;
const float FRICTION = 0.95;


///--------------------------------------------------------------------------
//---------------------------------------------------------------------------
///--------------------------------Structure---------------------------------
//---------------------------------------------------------------------------
///--------------------------------------------------------------------------
struct Complexe
{
    float x,y;
};


Complexe make_force(float x, float y)
{
    Complexe c;
    c.x=x;
    c.y=y;
    return c;
}

Complexe make_complex_exp(float r, float theta)
{
    Complexe c;
    c.x=r*cos(theta);
    c.y=r*sin(theta);
    return c;
}

Complexe operator+(Complexe a, Complexe b)
{
    Complexe c;
    c.x=a.x+b.x;
    c.y=a.y+b.y;
    return c;
}


Complexe operator-(Complexe a, Complexe b)
{
    Complexe c;
    c.x=a.x-b.x;
    c.y=a.y-b.y;
    return c;
}

Complexe operator*(Complexe a, Complexe b)
{
    Complexe c;
    c.x=a.x*b.x-a.y*b.y;
    c.y=a.x*b.y+a.y*b.x;
    return c;
}

Complexe operator*(Complexe a, float b)
{
    Complexe c;
    c.x=a.x*b;
    c.y=a.y*b;
    return c;
}

Complexe operator*(float b, Complexe a)
{
    Complexe c;
    c.x=a.x*b;
    c.y=a.y*b;
    return c;
}

typedef Complexe Vect2;
typedef Complexe Point;

float norme(Vect2 a)
{
    return sqrt((a.x*a.x)+(a.y*a.y));
}

float scalaire(Vect2 A, Vect2 B)
{
    return A.x*B.x+A.y*B.y;
}

Point make_pos(float x, float y)
{
    Point p;
    p.x=x;
    p.y=y;
    return p;
}

Vect2 make_vect(float x, float y)
{
    Vect2 v;
    v.x=x;
    v.y=y;
    return v;
}

struct Couleur
{
    int r,g,b;
};

Couleur make_couleur(int r, int g, int b)
{
    Couleur c;
    c.r = r;
    c.g = g;
    c.b = b;
    return c;
}


struct Particule
{
    Point p;
    Vect2 v;
    Vect2 f;
    int r;
    float m;
    Couleur c;
};

struct Ligne
{
    Point a;
    Point b;
    Vect2 vectAB;
    float m,p; ///Eq : y = m*x+p
    Couleur c;
};


struct Nouvelle_particule
{
    int x1,y1;
    int x2,y2;
    int etat; ///0: RIEN__1:En cours de création__2 : Terminé
    int r;
    float m;
    bool choc;
};

struct Nouvelle_ligne
{
    Point a;
    Point b;
    Vect2 vectAB;
    float m,p;
    int etat; ///0: RIEN__1:En cours de création__2 : Terminé
    bool choc;
    bool mouse; /// 0 : droite et 1 : gauche
};

struct World
{
    int n;
    int nbligne;

    Particule tab[MAX_PART];
    Ligne tabLg[MAX_LIGNE];

    bool affichage; ///Pour debug
    bool creation;  /// 0 : Ligne et 1 : Balle

    Nouvelle_particule new_part;
    Nouvelle_ligne new_line;
};

///--------------------------------------------------------------------------
//---------------------------------------------------------------------------
///------------------------COLLISON ENTRE BALLES-----------------------------
//---------------------------------------------------------------------------
///--------------------------------------------------------------------------

float distance_2pts(Point p1, Point p2)
{
    return sqrt((p1.x-p2.x)*(p1.x-p2.x)+((p1.y-p2.y)*(p1.y-p2.y)));
}

float collision(Point A, int ra, Point B, int rb)
{
    float d;
    d=distance_2pts(A,B);
    if(d<=ra+rb)
        return true;
    else
        return false;
}


void evite_chevauchement2(Point &A, int ra, Point &B, int rb)
{
    Vect2 unitaire;
    float distance_2centres;

    ///Distance entre les 2 centres
    distance_2centres = distance_2pts(A,B);

    ///Calcul du vecteur unitaire (Sens de B vers A)
    unitaire.x = (A.x-B.x) / distance_2centres;
    unitaire.y = (A.y-B.y) / distance_2centres;

    ///Reposition de la particule A grâce au vecteur unitaire
    A.x = B.x + (ra + rb + 1) * unitaire.x; ///Centrage des deux cercles --> Decalage (rayon A * rayon B) dans l'axe des x
    A.y = B.y + (ra + rb + 1) * unitaire.y; ///Centrage des deux cercles --> Decalage (rayon A * rayon B) dans l'axe des y
    ///+1 pour éviter que les deux cercles ne se confondent en 1 point
}

void evite_chevauchement(Particule &A, Particule &B,int i, int j, World &d)
{
    Vect2 unitaire;
    float distance_2centres;
    int tmp;

    ///Distance entre les 2 centres
    distance_2centres = distance_2pts(A.p,B.p);

    ///Calcul du vecteur unitaire (Sens de B vers A)
    unitaire.x = (A.p.x-B.p.x) / distance_2centres;
    unitaire.y = (A.p.y-B.p.y) / distance_2centres;

    ///Reposition de la particule A grâce au vecteur unitaire
    A.p.x = B.p.x + (A.r + B.r + 1) * unitaire.x; ///Centrage des deux cercles --> Decalage (rayon A * rayon B) dans l'axe des x
    A.p.y = B.p.y + (A.r + B.r + 1) * unitaire.y; ///Centrage des deux cercles --> Decalage (rayon A * rayon B) dans l'axe des y
    ///+1 pour éviter que les deux cercles ne se confondent en 1 point

    ///Ce test permet d'éviter un décalage hors de la fenêtre
    if(A.p.y <= A.r || A.p.y >= DIMH - A.r || A.p.x <= A.r || A.p.x >= DIMH - A.r)
    {
        ///Calcul du vecteur unitaire (Sens de A vers B)
        unitaire.x = (B.p.x-A.p.x) / distance_2centres;
        unitaire.y = (B.p.y-A.p.y) / distance_2centres;

        ///Reposition de la particule A grâce au vecteur unitaire
        B.p.x = A.p.x + (A.r + B.r + 1) * unitaire.x; ///Centrage des deux cercles --> Decalage (rayon A * rayon B) dans l'axe des x
        B.p.y = A.p.y + (A.r + B.r + 1) * unitaire.y; ///Centrage des deux cercles --> Decalage (rayon A * rayon B) dans l'axe des y
        ///+1 pour éviter que les deux cercles ne se confondent en 1 point
        tmp=j;
        j=i;
        i=j;

    }

    ///Boucle limitant le risque de déplacer une balle dans une autre
    int k;
    for(k=0;k<d.n;k++)
    {
        if(k != i && k !=j)
        {
            if(collision(d.tab[i].p,d.tab[i].r,d.tab[k].p,d.tab[k].r))
            {
                evite_chevauchement2(d.tab[k].p,d.tab[k].r,d.tab[i].p,d.tab[i].r);
            }
        }
    }

}

void nouvelle_vitesse(Particule &A, Particule &B)
{
    Vect2 n,t;
    Vect2 va_n_apres2, va_t_apres2, vb_n_apres2, vb_t_apres2;
    float va_n, va_t, vb_n, vb_t;
    float va_n_apres, vb_n_apres;
    float dd;

    ///Calcul vecteur unitaire n (normal à la tangente)
    dd = distance_2pts(A.p,B.p);
    n = make_vect( (B.p.x - A.p.x) / dd, (B.p.y - A.p.y) / dd );

    ///Calcul vecteur unitaire t (directeur de la tangente)
    t = make_vect(-n.y,n.x);

    ///n et t forment un plan orthogonal
    ///On projette les vecteurs va et vb dessus pour se ramener au choc elastique 1 dimension (celle de n)
    va_n = scalaire(n,A.v);
    va_t = scalaire(t,A.v);

    vb_n = scalaire(n,B.v);
    vb_t = scalaire(t,B.v);

    ///Calcul du choc elastique entre A et B découlant de la 2eme loi de newton dans le referenciel (n,t)
    ///Aucunes forces ne s'exercent sur l'axe t, donc va_t_apres = va_t donc on ne calcule rien.
    va_n_apres = ((va_n * (A.m-B.m)) + (2*B.m*vb_n)) / (A.m + B.m);
    vb_n_apres = ((vb_n * (B.m-A.m)) + (2*A.m*va_n)) / (A.m + B.m);

    ///Multiplication par les vecteurs unitaires pour retrouver des vecteurs
    va_n_apres2 = va_n_apres * n;
    va_t_apres2 = va_t * t;

    vb_n_apres2 = vb_n_apres * n;
    vb_t_apres2 = vb_t * t;

    ///Adition des deux coordonnées pour se retrouver dans le référenciel de base + ajout d'un coefficient de friction)
    A.v = (va_n_apres2 + va_t_apres2) * FRICTION;
    B.v = (vb_n_apres2 + vb_t_apres2) * FRICTION;

}

///--------------------------------------------------------------------------
//---------------------------------------------------------------------------
///------------------------COLLISON PAROIS ET MUR----------------------------
//---------------------------------------------------------------------------
///--------------------------------------------------------------------------

void collision_murs(Particule &p)
{
    if(p.p.y <= p.r)
    {
        p.p.y = (- p.p.y) + p.r * 2;
        p.v.y = (- p.v.y) * FRICTION;
    }
    else
    {
        if(p.p.y >= DIMH - p.r)
        {
            p.p.y = DIMH + DIMH - p.p.y - p.r *2;
            p.v.y = (- p.v.y) * FRICTION;
        }
    }
    if(p.p.x <= p.r)
    {
        p.p.x = (-p.p.x) + p.r * 2;
        p.v.x = (-p.v.x) * FRICTION;
    }
    else
    {
        if(p.p.x >= DIMW - p.r)
        {
            p.p.x = DIMW + DIMW - p.p.x - p.r *2;
            p.v.x = (-p.v.x) * FRICTION;
        }
    }



}

bool test_collision_parois(Particule &p, Ligne lg, Point &collision_probable)
{
    Vect2 projete, LgVersCentre, unitaire;

    ///Calcul du vecteur allant d'une extremite de la ligne au cercle + calcul vecteur unitaire pour pouvoir projeter
    LgVersCentre = make_vect(p.p.x-lg.a.x,p.p.y-lg.a.y);
    unitaire = lg.vectAB * (1/norme(lg.vectAB));

    ///Projection orthogonal du vecteur LgVersCentre sur la droite
    projete = scalaire(LgVersCentre,unitaire) * unitaire;


    ///Un point C appartient à un segment AB si et seulement 0 <= scalaire(AB,AC) <= scalaire (AB,AB)
    ///On remplace AC par "projete" est AB par "lg.vectAB"

    if(scalaire(lg.vectAB,projete) <= 0)
    {
        ///Le point est en dehors du segment, la collision la plus probable est donc l'extremite du segment
        collision_probable = lg.a;
    }
    else
    {
        if(scalaire(lg.vectAB,lg.vectAB) <= scalaire(lg.vectAB,projete))
        {
            ///Le point est en dehors du segment, la collision la plus probable est donc l'extremite du segment
            collision_probable = lg.b;
        }
        else
        {
            ///La collision probable est la projection du centre du cercle sur la ligne
            collision_probable = lg.a + projete;
        }
    }

    ///Si la distance enre le pts de collision probable et le cercle est inférieur au rayon --> Il y a collision
    return collision(collision_probable,0,p.p,p.r);
}

void nouvelle_vitesse_parois(Particule &p, Ligne lg, Point impact)
{
    Vect2 n, t, unitaire, va_t, va_n, va_n_apres;
    float dd;

    ///Calcul vecteur unitaire n (Vecteur normal à la tangente)
    dd = distance_2pts(p.p,impact);
    n = make_vect( (impact.x - p.p.x) / dd, (impact.y - p.p.y) / dd );

    ///Calcul vecteur unitaire t (directeur de la tangente)
    t = make_vect(-n.y,n.x);

    ///Projection du vecteur vitesse sur l'axe n et t
    va_n = scalaire(n,p.v) * n;
    va_t = scalaire(t,p.v) * t;

    ///On se ramène à une collision dans le plan (n,t) où t serait le vecteur directeur de la parois de collision
    ///Comme avec un mur classique, il suffit d'inverser la composante normal au mur (ici la composante n)
    va_n_apres = -1 * va_n;

    ///Adition des deux coordonnées pour se retrouver dans le référenciel de base
    p.v = va_n_apres + va_t; ///--> Pas de FRICTION pour ne pas arrêter les balles lorsqu'elles roulent sur une ligne
}


///--------------------------------------------------------------------------
//---------------------------------------------------------------------------
///------------------------FONCTION hors collision---------------------------
//---------------------------------------------------------------------------
///--------------------------------------------------------------------------
void reaction_eau(Particule &p)
{
    ///Formule volume sphere : 4*PI*r^3/3
    ///Formule volume sphère tronquée : PI * h^2/3 * (3r - h)
    Point tmp;
    Ligne lg;
    Vect2 stokes;

    ///Ligne de la surface de l'eau
    lg.a.x = DIMW/2;
    lg.a.y = SURFACE_EAU;
    lg.b.x = DIMW;
    lg.b.y = SURFACE_EAU;
    lg.c = make_couleur(255,0,123);
    lg.vectAB = make_vect(lg.b.x-lg.a.x,lg.b.y-lg.a.y);

    ///Test pour savoir si la balle est dans l'eau
    if(test_collision_parois(p,lg,tmp) || (p.p.y + p.r < SURFACE_EAU && p.p.x + p.r > DIMW/2))
    {
        float volume_balle_immerge, volume_hors_eau,h;
        int mvol_eau,facteur;
        Vect2 Pa;

        ///Volume balle
        volume_balle_immerge = 4 * M_PI * p.r * p.r * p.r /3;

        ///Test de submersion
        h = (p.p.y + p.r - lg.a.y);

        if( h > 0) ///Balle non submergé
        {
            volume_hors_eau = M_PI * h * h * (3 * p.r - h) /3;
            volume_balle_immerge = volume_balle_immerge - volume_hors_eau;
        }

        ///Masse volumique de l'eau (non pur)
        mvol_eau = 997;

        ///Calcul de la poussée d'archimède avec un facteur de division pour rendre le résultat réaliste (facteur trouve à tâtons)
        ///Formule Pa = - masse volumique * Volume du corps * accelaration champ de pesanteur
        facteur = 5000000;
        Pa = -mvol_eau * volume_balle_immerge/facteur * make_vect(0,-9.81);

        ///Calcul de la loi de Stokes (=Force de frottement d'un fluide sur une sphère)
        ///Formule F = -6 * PI * Rayon * Vitesse * Viscosité dynamique (0,0017 pour de l'eau proche de 0°C)
        stokes = -6*M_PI*p.r*p.v*0.0017;

        ///Actualisation des forces
        p.f = p.f + Pa + stokes;
    }
}

void supprimer(World &d, int indice)
{
    int i;

    if(d.creation)
    {
        for(i=indice;i<d.n-1;i++)
        {
            d.tab[i] = d.tab[i+1];
        }
        d.n += -1;
    }
    else
    {
        for(i=indice;i<d.nbligne-1;i++)
        {
            d.tabLg[i] = d.tabLg[i+1];
        }
        d.nbligne += -1;
    }

}

bool appartient_seg(Point C, Point A, Vect2 AB)
{
    Vect2 AC;
    float scal;

    ///Un point C appartient à un segment AB si et seulement 0 <= scalaire(AB,AC) <= scalaire (AB,AB)
    AC = make_vect(C.x-A.x,C.y-A.y);
    scal = scalaire(AB,AC);
    return ( 0 <= scal && scal <= scalaire(AB,AB) );

}

bool secante_bassin(float m, float p, Vect2 AB, Point A)
{
    Point C;
    Vect2 AC;
    float scal;

    if(AB.x == DIMW/2) ///Droites parallèles
    {
        return false;
    }
    else
    {
        ///Recherche du point d'intersection (C) entre la ligne et la parois du bassin
        C.x = DIMW/2;
        C.y = m * C.x + p;

        return ((C.y >=0 && C.y <= SURFACE_EAU) && appartient_seg(C,A,AB));
    }
}

bool secante(Ligne lg, Nouvelle_ligne new_lg)
{
    Point C;
    Vect2 AC, AD;
    float scal,scal2;
    bool retour;
    ///Il y a 4 cas possibles à étudier après avoir écarte que les 2 lignes sont strictement parallèle
    /// --> Ligne strictement parallèle
    /// --> Ligne sur la même droite
    /// --> Une des deux lignes est verticale
    /// --> Les 2 lignes sont obliques
    ///Methode : recherche du point d'intersection des droites puis on regarde s'il appartient aux segments


    if((lg.m == new_lg.m && lg.p != new_lg.p) || (lg.p == -1 && new_lg.p == -1 && lg.a.x != new_lg.a.x) ) ///Elimine le cas des lignes strictement parallèles
    {
        retour = false;
    }
    else
    {
        if(lg.m == new_lg.m) ///Ligne sur la même droite
        {
            ///On recherche si les extrémités du segment le plus petit appartiennent au segment le plus grand.
            if(norme(lg.vectAB) > norme(new_lg.vectAB))
            {
                retour = appartient_seg(new_lg.a,lg.a,lg.vectAB) || appartient_seg(new_lg.b,lg.a,lg.vectAB);
            }
            else
            {
                retour = appartient_seg(lg.a,new_lg.a,new_lg.vectAB) || appartient_seg(lg.b,new_lg.a,new_lg.vectAB);;
            }


        }
        else ///Une des deux ligne est verticale
        {
            if(lg.p == -1) ///Lg est verticale
            {
                C.x = lg.a.x;
                C.y = new_lg.m * C.x + new_lg.p;
                retour = appartient_seg(C, new_lg.a,new_lg.vectAB);
            }
            else
            {
                if(new_lg.p == -1) ///new_lg est verticale
                {
                    C.x = new_lg.a.x;
                    C.y = lg.m * C.x + lg.p;
                    retour = appartient_seg(C, lg.a,lg.vectAB);
                }
                else ///Lignes obliques
                {
                    ///Recherche du point d'intersection (C) entre les 2 lignes
                    C.x = (-lg.p + new_lg.p) / (lg.m - new_lg.m);
                    C.y = lg.m * C.x + lg.p;
                    retour = appartient_seg(C,lg.a,lg.vectAB) && appartient_seg(C,new_lg.a,new_lg.vectAB);
                }

            }
        }
    }
    return retour;
}

void addParticule(World &d)
{
    Particule p;
    if(d.n<MAX_PART)
    {
        p.p = make_pos(d.new_part.x1,d.new_part.y1);
        p.v = make_vect(d.new_part.x1-d.new_part.x2,d.new_part.y1-d.new_part.y2);
        p.f = make_force(0,0);
        p.r = d.new_part.r;
        p.m = d.new_part.m;
        p.c = make_couleur(irand(0,255),irand(0,255),irand(0,255));
        d.tab[d.n] = p;
        d.n = d.n + 1;
    }
}

void addLigne(World &d)
{
    ///Le vecteur nul n'est pas une ligne donc on l'élimine.
    if(d.nbligne<MAX_LIGNE && (d.new_line.vectAB.x != 0 || d.new_line.vectAB.y != 0) )
    {
        d.tabLg[d.nbligne].a = d.new_line.a;
        d.tabLg[d.nbligne].b = d.new_line.b;
        d.tabLg[d.nbligne].vectAB = d.new_line.vectAB;
        d.tabLg[d.nbligne].c = make_couleur(irand(0,255),irand(0,255),irand(0,255));
        d.tabLg[d.nbligne].m = d.new_line.m;
        d.tabLg[d.nbligne].p = d.new_line.p;
        d.nbligne++;
    }


}

void partAddForce(Particule &p)
{
    p.f = make_force(0,0);
    p.f = p.f + make_force(0, -p.m * 9.81);

    ///Poussée d'archimède
    reaction_eau(p);

}

void partUpdatePV(Particule &p)
{
    ///On divise par 15 pour rendre l'écart entre les rebonds de balle lourde et légère moins important
    p.v = p.v + (deltaT/15) * p.f;

    p.p = p.p + p.v*deltaT;
}

///--------------------------------------------------------------------------
//---------------------------------------------------------------------------
///---------------------FONCTIONS D'INITIALISATIONS--------------------------
///---------------------- + FONCTIONS PRINCIPALES ---------------------------
//---------------------------------------------------------------------------
///--------------------------------------------------------------------------
Particule PartInit(World &d)
{
    Particule part;
    int i,j;

    part.p = make_pos(irand(20,DIMW-50),irand(200,DIMH-50));
    part.r = irand(12,45);
    for(i=0;i<d.n;i++)
    {
        if(collision(part.p,part.r,d.tab[i].p,d.tab[i].r))
        {
            part.p = make_pos(irand(20,DIMW-50),irand(200,DIMH-50));
            part.r = irand(12,45);
            i=0;
        }
    }

    part.v = make_vect(irand(-100,100),irand(10,100));
    part.f = make_force(0,0);

    part.m = 4.0 * M_PI * part.r * part.r * part.r /FACTEUR_MASSE + 1.0 ; ///+1 pour éviter une masse de 0
    part.c = make_couleur(irand(0,255),irand(0,255),irand(0,255));
    return part;
}

Ligne LigneInit()
{
    Ligne lg;
    float m,p;

    ///On recommence le placement des lignes tant qu'elles ont un contact avec l'eau
    do
    {
        lg.a.x = irand(0,DIMW);
        lg.a.y = irand(0,DIMH);
        lg.b.x = irand(0,DIMW);
        lg.b.y = irand(0,DIMH);
        lg.vectAB = make_vect(lg.b.x-lg.a.x,lg.b.y-lg.a.y);
        if(lg.vectAB.x == 0) ///Cas spécial pour signaler une ligne horizontale
        {
            lg.m = 0;
            lg.p = -1;
        }
        else
        {
            lg.m = lg.vectAB.y/lg.vectAB.x;
            lg.p = lg.a.y - lg.m * lg.a.x;
        }
    }while((lg.a.x >= DIMW/2 && lg.a.y <= SURFACE_EAU) || (lg.b.x >= DIMW/2 && lg.b.y <= SURFACE_EAU) || (secante_bassin(lg.m,lg.p,lg.vectAB,lg.a)));

    lg.c = make_couleur(irand(0,255),irand(0,255),irand(0,255));

    return lg;
}

Ligne init_bassin()
{
    Ligne bassin;

    bassin.a.x=DIMW/2;
    bassin.a.y=0;
    bassin.b.x=DIMW/2;
    bassin.b.y=SURFACE_EAU;
    bassin.c = make_couleur(234,196,3);
    bassin.vectAB = make_vect(bassin.b.x-bassin.a.x,bassin.b.y-bassin.a.y);
    bassin.m = 0;
    bassin.p = -1;

    return bassin;
}

void init(World& d,int nbPart, int nbLigne)
{
    int i,j;

    d.nbligne = nbLigne + 1; ///+1 à cause de la ligne du bassin

    d.tabLg[0] = init_bassin();
    for(j=1;j<d.nbligne ;j++)
    {
        d.tabLg[j] = LigneInit();
    }

    d.n = 0;
    for(i=0;i<nbPart;i++)
    {
        d.tab[i]=PartInit(d);
        d.n += 1;
    }

    ///Initalisation des booléens
    d.affichage = false;
    d.creation = true;

    ///Initialisation de la nouvelle particule et de la nouvelle ligne
    d.new_part.etat = 0;
    d.new_part.r = 25;
    d.new_part.m = 4.0 * M_PI * d.new_part.r * d.new_part.r * d.new_part.r /FACTEUR_MASSE + 1.0;
    d.new_part.choc = false;

    d.new_line.etat = 0;
    d.new_line.choc = false;
    d.new_line.mouse = true;
}



void update(World &d)
{
    int i,j,k;
    int x,y;
    Point impact_probable;

    ///Contrôle des appuis sur touches du clavier
    if(isKeyPressed(SDLK_TAB))
        d.creation = !d.creation;
    if(isKeyPressed(SDLK_RETURN))
        d.affichage = !d.affichage;

    if(isKeyPressed(SDLK_UP) || isKeyPressed(SDLK_z))
    {
        if(d.new_part.r<=43)
            d.new_part.r += 2;
            d.new_part.m = 4.0 * M_PI * d.new_part.r * d.new_part.r * d.new_part.r /FACTEUR_MASSE + 1.0;
    }
    if(isKeyPressed(SDLK_DOWN) || isKeyPressed(SDLK_s))
    {
        if(d.new_part.r >= 14)
            d.new_part.r -= 2;
            d.new_part.m = 4.0 * M_PI * d.new_part.r * d.new_part.r * d.new_part.r /FACTEUR_MASSE + 1.0;
    }


    ///Gestion du mode BALLE
    if(d.creation)
    {
        ///Suppression des balles en cas de clique droit
        if(isMousePressed(SDL_BUTTON_RIGHT))
        {
           mousePos(x,y);
           for(i=0;i<d.n;i++)
           {
               if(collision(d.tab[i].p,d.tab[i].r, make_pos(x,y),0))
               {
                   supprimer(d,i);
               }
           }
        }
        ///Vérification des chocs de la balle en cours de création
        if(d.new_part.etat !=0)
        {
            d.new_part.choc = false;
            for(i=0;i<d.n;i++)
            {
                if(collision(make_pos(d.new_part.x1,d.new_part.y1),d.new_part.r,d.tab[i].p,d.tab[i].r))
                {
                    d.new_part.choc = true;
                    break;
                }
            }

            ///Ajout d'une balle s'il n'y a pas de choc
            if(d.new_part.etat == 2)
            {
                if(!d.new_part.choc)
                {
                    addParticule(d);
                }
                d.new_part.etat = 0;
            }
        }
    }
    else ///Gestion du mode LIGNE
    {
        ///Si une ligne est en cours de création
        if(d.new_line.etat !=0) ///Calcul des paramètres
        {
            d.new_line.choc = false;
            d.new_line.vectAB = make_vect(d.new_line.b.x-d.new_line.a.x,d.new_line.b.y-d.new_line.a.y);
            if(d.new_line.vectAB.x == 0) ///Cas d'une ligne verticale
            {
                d.new_line.m = 0;
                d.new_line.p = -1;
            }
            else
            {
                d.new_line.m = d.new_line.vectAB.y/d.new_line.vectAB.x;
                d.new_line.p = d.new_line.a.y - d.new_line.m * d.new_line.a.x;
            }


            if(d.new_line.mouse) ///Le dernier clique est gauche
            {
                ///Vérification des collisions avec le bassin d'eau
                for(i=0;i<d.n;i++)
                {
                    if((d.new_line.a.x >= DIMW/2 && d.new_line.a.y <= SURFACE_EAU) || (d.new_line.b.x >= DIMW/2 && d.new_line.b.y <= SURFACE_EAU)
                       || (secante_bassin(d.new_line.m,d.new_line.p,d.new_line.vectAB,d.new_line.a)))
                    {
                        d.new_line.choc = true;
                        break;
                    }
                }
                ///Ajout d'une ligne s'il n'y a pas de choc
                if(d.new_line.etat == 2)
                {
                    if(!d.new_line.choc)
                    {
                        addLigne(d);
                    }
                    d.new_line.etat = 0;
                }
            }
            else ///Le dernier clique est droite
            {
                if(d.new_line.etat == 2) ///Si clique relaché et que la ligne n'est pas un point --> suppression de la selection
                {
                    if(d.new_line.vectAB.x != 0 || d.new_line.vectAB.y != 0)
                    {
                       for(i=1;i<d.nbligne;i++)
                        {
                            if(secante(d.tabLg[i],d.new_line))
                            {
                                supprimer(d,i);
                                i=0;
                            }

                        }
                    }
                   d.new_line.etat = 0;
                }
            }
        }
    }


    ///Boucle pour déplacer les balles
    for(i=0;i<d.n;i++)
    {
        partAddForce(d.tab[i]);
        partUpdatePV(d.tab[i]);
        collision_murs(d.tab[i]);
    }

    ///Boucle pour gérer les collisions
    for(i=0;i<d.n;i++)
    {
        ///Avec les lignes
        for(k=0;k<d.nbligne;k++)
        {
            if(test_collision_parois(d.tab[i], d.tabLg[k],impact_probable))
            {
                ///ETAPE 1 : Reposition des cercles pour éviter le chevauchement
                evite_chevauchement2(d.tab[i].p,d.tab[i].r,impact_probable,0);

                ///ETAPE 2 : Le rebond
                nouvelle_vitesse_parois(d.tab[i],d.tabLg[k], impact_probable);
            }
        }
        ///Avec les autres balles
        for(j=i+1;j<d.n;j++)
        {
            if(collision(d.tab[i].p,d.tab[i].r,d.tab[j].p,d.tab[j].r))
            {
               ///ETAPE 1 : Reposition des cercles pour éviter le chevauchement
                evite_chevauchement(d.tab[j],d.tab[i],i,j,d);

                ///ETAPE 2 : Le rebond
                nouvelle_vitesse(d.tab[i],d.tab[j]);

            }
        }
    }
}

void draw(World& d)
{
    int i,j;
    int x,y;
    bool clique;
    clique = false;

    ///Dessine les balles
    for(i=0;i<d.n;i++)
    {
        color(d.tab[i].c.r,d.tab[i].c.g,d.tab[i].c.b);
        circle(d.tab[i].p.x,d.tab[i].p.y,d.tab[i].r);

        ///Affichage du DEBUG
        if(d.affichage)
        {
            ///Dessine les vecteurs vitesses
            color(255,0,0);
            line(d.tab[i].p.x,d.tab[i].p.y,d.tab[i].p.x+d.tab[i].v.x,d.tab[i].p.y+d.tab[i].v.y);

            ///Dessine les forces
            color(0,0,255);
            line(d.tab[i].p.x,d.tab[i].p.y,d.tab[i].p.x+d.tab[i].f.x,d.tab[i].p.y+d.tab[i].f.y);

            ///Masse de l'objet
            color(215, 206, 72);
            print(d.tab[i].p.x,d.tab[i].p.y,d.tab[i].m);
        }
    }

    ///Dessine les parois
    for(j=0;j<d.nbligne;j++)
    {
        color(d.tabLg[j].c.r,d.tabLg[j].c.g,d.tabLg[j].c.b);
        line(d.tabLg[j].a.x,d.tabLg[j].a.y,d.tabLg[j].b.x,d.tabLg[j].b.y);
    }

    ///DE L'EAU
    color(0, 255, 226,40);
    rectangleFill(DIMW/2,0,DIMW,SURFACE_EAU);

    ///Nombre de balle dans le monde
    color(87, 41, 5);
    rectangleFill(DIMW-170,DIMH-55,DIMW,DIMH);

    color(255, 173, 0);
    fontSize(15);
    print(DIMW-165,DIMH-30,"Nombre de balles : ");
    print(DIMW-30,DIMH-30,d.n);
    print(DIMW-165,DIMH-50,"Nombre de Lignes : ");
    print(DIMW-30,DIMH-50,d.nbligne - 1);

    ///Création d'objet (balle ou ligne)
    if(d.creation) ///GESTION DES BALLES
    {
        ///Remise à zéro du mode LIGNE dans le cas où l'utilisateur aurait changé de mode en cours de création
        if(d.new_line.etat !=0)
            d.new_line.etat = 0;

        color(255, 69, 0);
        fontSize(12);
        print(10,DIMH-20,"Taper sur TABULATION pour changer de mode : ");
        color(240, 255, 0);
        print(270,DIMH-20, "BALLE");
        color(255, 173, 0);
        print(10,DIMH-40,"Ajouter une balle avec le clique gauche");
        print(10,DIMH-60,"Supprimer une balle avec le clique droit");
        if(!(d.n >= MAX_PART))
        {
            if (d.new_part.etat == 1) ///affichage balle en cours de création
            {
                color( 163, 166, 149 );
                ///Vérifie les collisions
                if(d.new_part.choc)
                {
                    color(255,0,0);
                    print(d.new_part.x1 - d.new_part.r - 5, d.new_part.y1 + d.new_part.r + 5,"COLLISION !");
                }
                else
                {
                    print(d.new_part.x1 - d.new_part.r - 110, d.new_part.y1 + d.new_part.r + 5,
                      "Modifier la taille du rayon avec les flèches HAUT et BAS ");
                }
                line(d.new_part.x1,d.new_part.y1,d.new_part.x2,d.new_part.y2);
                circle(d.new_part.x1,d.new_part.y1,d.new_part.r);
                print(d.new_part.x1,d.new_part.y1,d.new_part.r);

            }
            if(isMousePressed(SDL_BUTTON_LEFT)) ///Actualisation de la balle en création
            {
                int x,y;
                mousePos(x,y);
                if (d.new_part.etat!=1)
                {
                    d.new_part.etat = 1;
                    d.new_part.x1 = x;
                    d.new_part.y1 = y;
                }
                else
                {
                    d.new_part.x2 = x;
                    d.new_part.y2 = y;
                }
            }
            else
            {
                if (d.new_part.etat==1)
                {
                    d.new_part.etat=2; ///Etat : 2 --> Fin de création
                }

            }
        }
        else
        {
            ///Message d'alerte si l'on atteint le nombre maximum de balles
            color(255, 0, 0);
            fontSize(15);
            print(10,DIMH-80,"Nombre de balles maximum atteint ! ");
        }
    }
    else ///GESTION DES LIGNES
    {
        ///Remise à zéro du mode balle dans le cas où l'utilisateur aurait changé de mode en cours de création
        if(d.new_part.etat !=0)
            d.new_part.etat = 0;

        color(255, 69, 0);
        fontSize(12);
        print(10,DIMH-20,"Taper sur TABULATION pour changer de mode : ");
        color(240, 255, 0);
        print(270,DIMH-20, "LIGNE");
        color(255, 173, 0);
        print(10,DIMH-40,"Dessiner une ligne avec le clique gauche");
        print(10,DIMH-60,"Supprimer une ligne en faisant glisser la souris avec le clique droit");

        ///Affichage ligne en cours de création
        if (d.new_line.etat == 1)
        {
            color( 163, 166, 149 );
            ///Vérifie les collisions
            if(d.new_line.choc && d.new_line.mouse)
            {
                color(255,0,0);
                print(d.new_line.a.x - 50, d.new_line.a.y + 5,"COLLISION BASSIN");
            }
            else
            {
                if(d.new_line.mouse)
                {
                    print(d.new_line.a.x - 75, d.new_line.a.y + 5,"Tracer la ligne avec votre souris");
                }
                else
                {
                    print(d.new_line.a.x - 85, d.new_line.a.y + 5,"Glisser votre souris vers les lignes à supprimer");
                }

            }
            line(d.new_line.a.x, d.new_line.a.y,d.new_line.b.x, d.new_line.b.y);
        }

        ///Si l'on peut encore créer des lignes
        if(!(d.nbligne >= MAX_LIGNE))
        {
            if(isMousePressed(SDL_BUTTON_LEFT)) ///Gestion clique gauche
            {
                int x,y;
                mousePos(x, y);
                d.new_line.mouse = true;
                if (d.new_line.etat!=1)
                {
                    d.new_line.etat = 1;
                    d.new_line.a.x = x;
                    d.new_line.a.y = y;
                    d.new_line.b.x = x;
                    d.new_line.b.y = y;
                }
                else
                {
                    d.new_line.b.x = x;
                    d.new_line.b.y = y;
                }
                clique = true;
            }
        }
        else
        {
            ///Message d'alerte si l'on atteint le nombre maximum de lignes
            color(255, 0, 0);
            fontSize(15);
            print(10,DIMH-80,"Nombre de lignes maximum atteint ! ");
        }

        ///Gestion clique droit
        if(isMousePressed(SDL_BUTTON_RIGHT))
        {
            int x,y;
            mousePos(x, y);
            d.new_line.mouse = false;
            if (d.new_line.etat!=1)
            {
                d.new_line.etat = 1;
                d.new_line.a.x = x;
                d.new_line.a.y = y;
                d.new_line.b.x = x;
                d.new_line.b.y = y;
            }
            else
            {
                d.new_line.b.x = x;
                d.new_line.b.y = y;
            }
        }
        else ///Gestion pas de clique
        {
            if(d.new_line.etat == 1 && !clique)
            {
                d.new_line.etat = 2; ///Etat : 2 --> Fin de création
            }
        }
    }
}

int main(int , char** )
{
    World dat;

    bool stop=false;
	winInit("MyProg", DIMW, DIMH);
    backgroundColor( 22, 22, 22 );
    srand(time(NULL));
    init(dat,10,2);

    Menu menu;
    menu_add(menu, "Init");
    menu_add(menu, "Run");
    menu_setSelect(menu,1);

	while( !stop )
    {
        winClear();

        update(dat);
        draw(dat);

        if(menu_select(menu)==0)
        {
            init(dat,10,2);
            menu_setSelect(menu,1);
        }

        menu_draw(menu);
        stop = winDisplay();
    }
    winQuit();
	return 0;
}



