#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
from pyssa.gui.data_structures import settings
from gui.data_structures import safeguard


class SettingsSafeguard(safeguard.Safeguard):

    def __init__(self, settings_obj):
        # var which contains a settings object
        super().__init__()
        self.settings: settings.Settings = settings_obj

