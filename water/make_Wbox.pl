#! /usr/bin/perl

use strict;

my %Atoms;
my $resn = 1;
my $atomn = 1;
my $size_x = $ARGV[0];
my $size_y = $ARGV[1];
my $size_z = $ARGV[2];
my $x;
my $y;
my $z;
my $nx;
my $ny;
my $nz;
my $x_curr = 1;
my $y_curr = 1;
my $z_curr = 1;

open (box, ">box.pdb");

# determine number of water molecules #
# density of water is g/cm^3#
# Mass of water = 2.8971^-23 g/molecule #
# We have a volume of (in Angstrom m^3): #
my $volume = ($size_x*1E-8) * ($size_y*1E-8) * ($size_z*1E-8) ;

#Therefore we need molecules:
my $mass = $volume * 1;
my $nmolecules = $mass / 2.8971E-23 ;
#print "$nmolecules\n";
$volume = ($size_x) * ($size_y) * ($size_z);
$nx = $size_x*(($nmolecules/$volume)**(1/3));
$ny = $size_y*(($nmolecules/$volume)**(1/3));
$nz = $size_z*(($nmolecules/$volume)**(1/3));
#print "$nx, $ny, $nz\n";

my $one_side_x = $size_x/$nx;
my $one_side_y = $size_y/$ny;
my $one_side_z = $size_z/$nz;
#print "Space for $nmolecules molecule: side length: $one_side_x, $one_side_y,$one_side_z\n";

$x = -1*($size_x/2);
while($x_curr <= $nx){
  $y_curr = 1;
  $y = -1*($size_y/2);
  while($y_curr <= $ny){
    $z_curr = 1;
    $z = -1*($size_z/2);
    while($z_curr <= $nz){
      $Atoms{'O'}{$resn}{$atomn}{'x'} = $x;
      $Atoms{'O'}{$resn}{$atomn}{'y'} = $y;
      $Atoms{'O'}{$resn}{$atomn}{'z'} = $z;
      $Atoms{'O'}{$resn}{$atomn}{type} = "ATOM";
      $Atoms{'O'}{$resn}{$atomn}{name} = "OH2";
      $Atoms{'O'}{$resn}{$atomn}{resname} = "TIP3";
      $Atoms{'O'}{$resn}{$atomn}{occup} = 1.00;
      $Atoms{'O'}{$resn}{$atomn}{temp} = 0.00;

      $Atoms{'H1'}{$resn}{$atomn+1}{'x'} = $x -0.14;
      $Atoms{'H1'}{$resn}{$atomn+1}{'y'} = $y -0.76;
      $Atoms{'H1'}{$resn}{$atomn+1}{'z'} = $z -0.56;
      $Atoms{'H1'}{$resn}{$atomn+1}{type} = "ATOM";
      $Atoms{'H1'}{$resn}{$atomn+1}{name} = "H1";
      $Atoms{'H1'}{$resn}{$atomn+1}{resname} = "TIP3";
      $Atoms{'H1'}{$resn}{$atomn+1}{occup} = 1.00;
      $Atoms{'H1'}{$resn}{$atomn+1}{temp} = 0.00;

      $Atoms{'H2'}{$resn}{$atomn+2}{'x'} = $x -0.295;
      $Atoms{'H2'}{$resn}{$atomn+2}{'y'} = $y +1.08;
      $Atoms{'H2'}{$resn}{$atomn+2}{'z'} = $z -0.52;
      $Atoms{'H2'}{$resn}{$atomn+2}{type} = "ATOM";
      $Atoms{'H2'}{$resn}{$atomn+2}{name} = "H2";
      $Atoms{'H2'}{$resn}{$atomn+2}{resname} = "TIP3";
      $Atoms{'H2'}{$resn}{$atomn+2}{occup} = 1.00;
      $Atoms{'H2'}{$resn}{$atomn+2}{temp} = 0.00;








      
      printf box ("%-6s%5d% 4s% 6s%2s%4d%11.3f%8.3f%8.3f%6.2f%6.2f\n", $Atoms{'O'}{$resn}{$atomn}{type}, $atomn, $Atoms{'O'}{$resn}{$atomn}{name},$Atoms{'O'}{$resn}{$atomn}{resname},"",$resn, $Atoms{'O'}{$resn}{$atomn}{'x'}, $Atoms{'O'}{$resn}{$atomn}{'y'}, $Atoms{'O'}{$resn}{$atomn}{'z'},$Atoms{'O'}{$resn}{$atomn}{occup},$Atoms{'O'}{$resn}{$atomn}{temp});
      $atomn++;
      printf box ("%-6s%5d% 4s% 6s%2s%4d%11.3f%8.3f%8.3f%6.2f%6.2f\n", $Atoms{'H1'}{$resn}{$atomn}{type}, $atomn, $Atoms{'H1'}{$resn}{$atomn}{name},$Atoms{'H1'}{$resn}{$atomn}{resname},"",$resn, $Atoms{'H1'}{$resn}{$atomn}{'x'}, $Atoms{'H1'}{$resn}{$atomn}{'y'}, $Atoms{'H1'}{$resn}{$atomn}{'z'},$Atoms{'H1'}{$resn}{$atomn}{occup},$Atoms{'H1'}{$resn}{$atomn}{temp});
      $atomn++;
      printf box ("%-6s%5d% 4s% 6s%2s%4d%11.3f%8.3f%8.3f%6.2f%6.2f\n", $Atoms{'H2'}{$resn}{$atomn}{type}, $atomn, $Atoms{'H2'}{$resn}{$atomn}{name},$Atoms{'H2'}{$resn}{$atomn}{resname},"",$resn, $Atoms{'H2'}{$resn}{$atomn}{'x'}, $Atoms{'H2'}{$resn}{$atomn}{'y'}, $Atoms{'H2'}{$resn}{$atomn}{'z'},$Atoms{'H2'}{$resn}{$atomn}{occup},$Atoms{'H2'}{$resn}{$atomn}{temp});


      $z = $z + $one_side_z;
      $z_curr++;
      $resn++;
      $atomn++;
    };
    $y = $y + $one_side_y;
    $y_curr++;
  }; 
  $x = $x + $one_side_x;
  $x_curr++;
};
print box "\n";
$atomn = ($atomn) / 3;
printf ("placed %d in box instead of %d molecules. Error : %.3f \%\n", $atomn, $nmolecules, (1- $atomn/$nmolecules)*100);
