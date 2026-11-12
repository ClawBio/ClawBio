// Progressive enhancement: all skill content is already present in the HTML.
const search = document.getElementById('search');
const rows = [...document.querySelectorAll('#results article')];
const searchable = rows.map(row => (row.textContent + ' ' + (row.dataset.searchKeywords || '')).toLowerCase());
function update() {
  const words = search.value.toLowerCase().trim().split(/\s+/).filter(Boolean);
  let count = 0;
  rows.forEach((row, index) => {
    row.hidden = !words.every(word => searchable[index].includes(word));
    if (!row.hidden) count++;
  });
  document.getElementById('count').textContent = `${count} of ${rows.length} skills`;
  document.getElementById('empty').hidden = count !== 0;
}
search.addEventListener('input', update);
update();
