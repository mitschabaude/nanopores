import fs from 'fs';
import {exec} from 'child_process';

let colors = {
  pnps: 'ff8000',
  pnpsDark: 'a25100',
  protein: '1f77b4',
  proteinDark: '165480',
  receptor: 'ff99bb',
  experiment: 'ff0000',
  simulations: '00cc00',
};

let rotateRed = blue => blue.slice(4, 6) + blue.slice(0, 4);
let rotateGreen = blue => blue.slice(2, 6) + blue.slice(0, 2);

let light = '788abd';
let medium = '3c59bd';
let intense = '0224bd';

let complement = 'c2b342'; // rather ugly yellowish
let gold = 'c7a708';
let orange = 'cc9c43';

let newColors = {
  experiment: rotateRed(medium),
  simulations: medium,

  protein: medium,
  proteinDark: medium,

  pnps: gold,
  pnpsDark: gold,

  receptor: gold,
};

let svg = fs.readFileSync('./tmp/figures/fig1/first_figure_new.svg', {
  encoding: 'utf-8',
});

for (let c in colors) {
  if (c in newColors) {
    svg = svg.replace(new RegExp(colors[c], 'g'), newColors[c]);
  }
}

fs.writeFileSync('./tmp/figures/fig1/first_figure.svg', svg);
exec(
  'inkscape tmp/figures/fig1/first_figure.svg --export-pdf=./tmp/figures/fig1.pdf'
);
