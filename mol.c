
#include "mol.h"
#include<string.h>
#include <stdlib.h>




void atomset(atom *a, char element[3], double *x, double *y, double *z)
{
   strncpy(a->element, element, 3);
    a->x = *x;
    a->y = *y;
    a->z = *z;
}

void atomget(atom *a, char element[3], double *x, double *y, double *z)
{
    strcpy(element, a->element);
    *x = a->x;
    *y = a->y;
    *z = a->z;
}
// Function to set the values of a bond
void bondset( bond *bond, atom *a1, atom *a2, unsigned char epairs )
{
    bond->a1 = a1;
    bond->a2 = a2;
    bond->epairs = epairs;
}
// Function to get the values of a bond
void bondget(bond *b, atom **a1, atom **a2, unsigned char *epairs)
{
    *a1 = b->a1;
    *a2 = b->a2;
    *epairs = b->epairs;
}
molecule *molmalloc(unsigned short atom_max, unsigned short bond_max)
{
    molecule *m = malloc(sizeof(molecule));
    if (m == NULL)
    {
        return NULL;
    }
    m->atom_max = atom_max;
    m->atom_no = 0;
    m->atoms = malloc(sizeof(atom) * atom_max);
    if (m->atoms == NULL)
    {
        free(m);
        return NULL;
    }
    m->atom_ptrs = malloc(sizeof(atom *) * atom_max);
    if (m->atom_ptrs == NULL)
    {
        free(m->atoms);
        free(m);
        return NULL;
    }
    m->bond_max = bond_max;
    m->bond_no = 0;
    m->bonds = malloc(sizeof(bond) * bond_max);
    if (m->bonds == NULL)
    {
        free(m->atoms);
        free(m->atom_ptrs);
        free(m);
        return NULL;
    }
    m->bond_ptrs = malloc(sizeof(bond *) * bond_max);
    if (m->bond_ptrs == NULL)
    {
        free(m->atoms);
        free(m->atom_ptrs);
        free(m->bonds);
        free(m);
        return NULL;
    }
    return m;
}

molecule *molcopy(molecule *src)
{
    molecule *mol = molmalloc(src->atom_max, src->bond_max);

    mol->atom_no = src->atom_no;
    for(int i = 0; i < mol->atom_no; i++)
    {
        atomset(&mol->atoms[i], src->atoms[i].element, &src->atoms[i].x, &src->atoms[i].y, &src->atoms[i].z);
        mol->atom_ptrs[i] = &mol->atoms[i];
    }

    mol->bond_no = src->bond_no;
    for(int i = 0; i < mol->bond_no; i++)
    {
      bondset(&mol->bonds[i], mol->atom_ptrs[src->bonds[i].a1 - src->atoms], mol->atom_ptrs[src->bonds[i].a2 - src->atoms], src->bonds[i].epairs);
    }
    return mol;
}


void molfree(molecule *ptr)
{
    free(ptr->atoms);
    free(ptr->atom_ptrs);
    free(ptr->bonds);
    free(ptr->bond_ptrs);
    free(ptr);
}
void molappend_atom(molecule *molecule, atom *atom)
{
    if (molecule->atom_no == molecule->atom_max)
    {
        if (molecule->atom_max == 0)
        {
            molecule->atom_max = 1;
        }
        else
        {
            molecule->atom_max *= 2;
        }
        molecule->atoms = (struct atom *)realloc(molecule->atoms, sizeof(struct atom) * molecule->atom_max);
        molecule->atom_ptrs = (struct atom **)realloc(molecule->atom_ptrs, sizeof(struct atom *) * molecule->atom_max);
    }
    molecule->atoms[molecule->atom_no] = *atom;
    molecule->atom_ptrs[molecule->atom_no] = &molecule->atoms[molecule->atom_no];
    molecule->atom_no++;
}
void molappend_bond(molecule *molecule, bond *bond)
{
    if (molecule->bond_no == molecule->bond_max)
    {
        if (molecule->bond_max == 0)
            molecule->bond_max = 1;
        else
            molecule->bond_max *= 2;
        molecule->bonds = (struct bond *)realloc(molecule->bonds, sizeof(bond) * molecule->bond_max);
        molecule->bond_ptrs = (struct bond **)realloc(molecule->bond_ptrs, sizeof(struct bond *) * molecule->bond_max);
    }
    molecule->bonds[molecule->bond_no] = *bond;
    molecule->bond_ptrs[molecule->bond_no] = &molecule->bonds[molecule->bond_no];
    molecule->bond_no++;
}

// Compare function for atoms
    int atom_cmp(const void *a, const void *b)
    {
        atom *aa = *(atom **)a;
        atom *bb = *(atom **)b;
        return (aa->z > bb->z) - (aa->z < bb->z);
    }

    // Compare function for bonds
    int bond_cmp(const void *a, const void *b)
    {
        bond *ab = *(bond **)a;
        bond *bb = *(bond **)b;
        float az = (ab->a1->z + ab->a2->z) / 2;
        float bz = (bb->a1->z + bb->a2->z) / 2;
        return (az > bz) - (az < bz);
    }


void molsort(molecule *mol)
{
    qsort(mol->atom_ptrs, mol->atom_no, sizeof(atom *), atom_cmp);
    qsort(mol->bond_ptrs, mol->bond_no, sizeof(bond *), bond_cmp);
}
